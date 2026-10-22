"""Pillar D: estimate the trading policy rather than the forecast.

A linear policy holds, in each asset, a position that is a lag curve applied
to the asset's own past scaled returns (with the volatility-state surface),
divided by volatility. Its daily profit is ``beta' p_t`` with
``p_t = sum_i X_{i,t}' r_{i,t+1}``, the design rows of all assets weighted by
their next scaled return. Mean-variance utility with a quadratic cost of
trading is quadratic in ``beta``, so ``exp(utility)`` is a Gaussian
pseudo-likelihood and the whole collapsed machinery applies: it is the
regression of a vector of ones on the rows ``p_t`` (Britten-Jones, 1999),
plus pseudo-rows ``sqrt(c) (X_t - X_{t-1})`` with response zero for the
cost, under the adaptive smoothness prior. Unlike the forecast-based rule
the objective sees the portfolio, so redundant signals across correlated
assets are discounted automatically. The policy is compared walk-forward
with the forecast-based rule on identical rows.

Run: python examples/policy_curves.py [--final] [--cost 0.0005] [--risk 1.0]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel import HOLDOUT_START, load_wide_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows  # noqa: E402
from factor_curves import DROP  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import FastFitter, subperiod_sharpe, walk_forward  # noqa: E402
from badlm.design import LagComponent, SuffStats, build_multi_design  # noqa: E402
from badlm.model import build_modifier_bases  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")


def policy_walk_forward(model, df, modifier, start, step, cost_weight, risk):
    """Walk-forward estimate of the policy curve; returns per-row positions
    in scaled units (to be divided by volatility) and the curves."""
    model.standardize = False
    x = df["x"].to_numpy()[:, None]
    mod_bases = build_modifier_bases(modifier, x, model.modifier_k, 1, standardize=False) if model.modifier_k else [None]
    model.components_ = [LagComponent(model.structs_[0], mod_bases[0])]
    model.C_ = None
    data, _ = build_multi_design(x, df["y"].to_numpy(), model.components_, intercept=False, standardize=False,
                                 groups=df["asset"].to_numpy(), modifier=modifier)
    X, y = data.X, data.y
    t_row = df["time"].to_numpy()[data.keep]
    a_row = df["asset"].to_numpy()[data.keep]
    # profit rows p_t = sum_i X_it y_it, one per day
    days, inv = np.unique(t_row, return_inverse=True)
    P = np.zeros((days.size, X.shape[1]))
    np.add.at(P, inv, X * y[:, None])
    # cost rows: change of each asset's design row from one day to the next
    order = np.lexsort((t_row, a_row))
    Xo, ao, to = X[order], a_row[order], t_row[order]
    same = ao[1:] == ao[:-1]
    dX = (Xo[1:] - Xo[:-1])[same]
    dt = to[1:][same]
    fitter = None
    pos = np.full(y.size, np.nan)
    curves, T = [], start
    while T < t_row.max():
        tr_days = days <= T
        tr_cost = dt <= T
        test = (t_row > T) & (t_row <= T + step)
        if tr_days.sum() > 200 and test.any():
            rows_X = np.vstack([np.sqrt(risk) * P[tr_days], np.sqrt(cost_weight) * dX[tr_cost]])
            rows_y = np.r_[np.ones(tr_days.sum()) / np.sqrt(risk), np.zeros(tr_cost.sum())]
            stats = SuffStats.from_data(rows_X, rows_y)
            if fitter is None:
                fitter = FastFitter(model, stats, 0)
            res = fitter.fit(stats)
            pos[test] = X[test] @ res["mean"]
            curves.append(model.components_[0].curve(res["mean"], weights=None) if model.modifier_k is None
                          else model.components_[0].curve(res["mean"], z=model.components_[0].modifier.grid(5)))
        T += step
    have = np.isfinite(pos)
    return {"row": data.keep[have], "pos": pos[have], "curves": curves}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--cost", type=float, default=0.0005)
    ap.add_argument("--risk", type=float, default=1.0)
    args = ap.parse_args()
    panel = load_wide_panel()
    close, kind = panel["close"], panel["kind"]
    close = close.drop(columns=[c for c in DROP if c in close])
    if not args.final:
        close = close.loc[:HOLDOUT_START - pd.Timedelta(days=1)]
    df = stacked_rows(close, kinds=kind, use=("equity_us", "equity_uk", "index", "fx", "commodity"))
    med = df.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    df["volstate"] = np.log(df["vol"] / med).fillna(0.0).clip(-2, 2)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    breaks = [int(close.index.get_indexer([pd.Timestamp(y)], method="bfill")[0]) for y in ("2004-01-01", "2012-01-01")]
    pri = Priors(tau_level_scale=0.001, temper=10)
    print(f"{df['asset'].nunique()} assets, {len(df)} rows", flush=True)

    t0 = time.time()
    wf = walk_forward(BADLM(nlag=250, k=25, modifier_k=5, priors=pri), df["x"].to_numpy(), df["y"].to_numpy(),
                      df["time"].to_numpy(), groups=df["asset"].to_numpy(), start=start, step=args.step,
                      modifier=df["volstate"].to_numpy())
    rows = df.iloc[wf.row].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf.mean, wf.var
    results = evaluate("forecast rule (surface)", (rows["pred"] / rows["pvar"]).to_numpy(), rows,
                       rules=(("aim", {}), ("gp", {"speed": 0.05})))
    print(f"forecast rule: IC {ic(rows['pred'], rows['y']):.4f} ({time.time() - t0:.0f}s)", flush=True)
    out = {"forecast": [{k: v for k, v in r.items() if k not in ("pnl", "traded", "time")} for r in results]}
    for cw, label in ((0.0, "policy, no cost term"), (args.cost * 2000.0, f"policy, cost term {args.cost}")):
        t0 = time.time()
        pw = policy_walk_forward(BADLM(nlag=250, k=25, modifier_k=5, priors=pri), df, df["volstate"].to_numpy(),
                                 start, args.step, cost_weight=cw, risk=args.risk)
        prow = df.iloc[pw["row"]].reset_index(drop=True)
        res = evaluate(label, pw["pos"], prow, rules=(("aim", {}), ("gp", {"speed": 0.05})))
        slow = res[1]
        sub = subperiod_sharpe(slow["pnl"] - 0.0005 * slow["traded"], slow["time"], breaks)
        c = pw["curves"][-1]
        print(f"{label}: positions IC with outcome {np.corrcoef(pw['pos'], prow['y'])[0, 1]:.4f}; sub-periods "
              f"{[round(s[2], 2) for s in sub]}; final curve at mid state: lag1 {c[2][0]:+.4f} lag5 {c[2][4]:+.4f} "
              f"lag21 {c[2][20]:+.4f} sum {c[2].sum():+.3f} ({time.time() - t0:.0f}s)", flush=True)
        results += res
        out[label] = [{k: v for k, v in r.items() if k not in ("pnl", "traded", "time")} for r in res]
    print_table(results, "Forecast rule against estimated policy (development period)")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"policy_curves{'_final' if args.final else ''}.json"), "w") as f:
        json.dump({"args": vars(args), **out}, f, indent=1)
    print("written")
