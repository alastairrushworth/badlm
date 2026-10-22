"""Experiment 1: a pooled cross-asset lag curve as a trading signal.

Daily log returns of every series in the panel are divided by their ex-ante
volatility, stacked with one group per asset, and a single adaptive lag curve
of the next day's scaled return on the previous 250 days' scaled returns is
estimated. The curve is re-estimated at a fixed cadence on an expanding window
and each day's prediction uses only past data (``badlm.backtest.walk_forward``).
Predictions become risk-scaled positions, traded fully or slowly towards the
aim, and are evaluated gross and net of proportional costs against time-series
momentum and reversal rules on the same rows.

Run: python examples/trading_pooled.py [--universe all|fx|equity|commodity]
     [--k 84] [--scale 1.0] [--temper 1.0] [--step 63] [--start 1996] [--tag name]
"""
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _markets import load_panel  # noqa: E402
from _trading import benchmark_features, evaluate, ic, print_table, stacked_rows, strip, yearly_sharpe  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import deflated_sharpe, walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
resdir = os.path.join(here, "..", "paper", "results")


def parse():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="all")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--nlag", type=int, default=250)
    ap.add_argument("--k", type=int, default=84)
    ap.add_argument("--scale", type=float, default=1.0, help="tau_level_scale of the prior")
    ap.add_argument("--temper", type=float, default=1.0)
    ap.add_argument("--lag-scale", default="linear")
    ap.add_argument("--tag", default=None)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse()
    tag = args.tag or f"pooled_{args.universe}"
    panel = load_panel()
    close, kinds = panel["close"], panel["kind"]
    use = ("equity", "fx", "commodity") if args.universe == "all" else (args.universe,)
    df = stacked_rows(close, kinds=kinds, use=use)
    feats = benchmark_features(df)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    print(f"{df['asset'].nunique()} assets, {len(df)} rows; first refit at {close.index[start].date()}, "
          f"step {args.step} days, nlag {args.nlag}, k {args.k}, scale {args.scale}, temper {args.temper}", flush=True)

    model = BADLM(nlag=args.nlag, k=args.k, lag_scale=args.lag_scale,
                  priors=Priors(tau_level_scale=args.scale, temper=args.temper))
    t0 = time.time()
    wf = walk_forward(model, df["x"].to_numpy(), df["y"].to_numpy(), df["time"].to_numpy(),
                      groups=df["asset"].to_numpy(), start=start, step=args.step)
    print(f"walk-forward: {len(wf.refit_times)} refits in {time.time() - t0:.0f}s; "
          f"{wf.mean.size} out-of-sample predictions", flush=True)
    rows = df.iloc[wf.row].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf.mean, wf.var
    fb = feats.iloc[wf.row].reset_index(drop=True)
    oos_r2 = 1.0 - np.sum((rows["y"] - rows["pred"]) ** 2) / np.sum(rows["y"] ** 2)
    print(f"out-of-sample R^2 {100*oos_r2:.3f}%, IC {ic(rows['pred'], rows['y']):.4f}; benchmarks: "
          + ", ".join(f"{c} {ic(fb[c], rows['y']):.4f}" for c in ("sum5", "sum21", "sum63", "sum252"))
          + f"; period {close.index[rows['time'].min()].date()} to {close.index[rows['time'].max()].date()}")

    results = evaluate("badlm curve", (rows["pred"] / rows["pvar"]).to_numpy(), rows)
    combo = fb["sum21"] / np.sqrt(21) + fb["sum252"] / np.sqrt(252)
    for name, aim in (("TSMOM sign(252d)", np.sign(fb["sum252"])), ("TSMOM linear(252d)", fb["sum252"] / 252),
                      ("sum21 + sum252", combo), ("reversal -sum(5d)", -fb["sum5"] / 5)):
        results += evaluate(name, aim.to_numpy(), rows, rules=(("aim", {}), ("gp", {"speed": 0.1})))
    print_table(results, f"{args.universe}: out of sample, {close.index[rows['time'].min()].year}-{close.index[rows['time'].max()].year}")

    by_kind = {}
    if args.universe == "all":
        for kind in use:
            sel = (rows["asset"].map(kinds) == kind).to_numpy()
            sub = rows[sel].reset_index(drop=True)
            res = evaluate("badlm curve", (sub["pred"] / sub["pvar"]).to_numpy(), sub, rules=(("aim", {}), ("gp", {"speed": 0.1})))
            res += evaluate("TSMOM linear(252d)", (fb.loc[sel, "sum252"] / 252).to_numpy(), sub, rules=(("gp", {"speed": 0.1}),))
            print_table(res, f"{kind} only")
            by_kind[kind] = strip(res)

    n_trials = 3
    best = max(results[:3], key=lambda r: r["net_5bp"])
    pnl_net = best["pnl"] - 0.0005 * best["traded"]
    dsr, sr0 = deflated_sharpe(pnl_net, n_trials)
    print(f"\nbest badlm rule at 5 bp: {best['rule']}, net SR {best['net_5bp']:.2f}; deflated Sharpe probability "
          f"{dsr:.3f} (expected maximum of {n_trials} unskilled rules: {sr0*np.sqrt(252):.2f} annualised)")
    ys = yearly_sharpe(best, close.index)
    print("yearly net Sharpe (5 bp):", " ".join(f"{y}:{s:.1f}" for y, s in ys.items()))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    ax = axes[0]
    curves = np.array([c[0] for c in wf.curves])
    lags = np.arange(1, args.nlag + 2)
    for i, c in enumerate(curves):
        ax.plot(lags, c, color=plt.cm.viridis(i / max(len(curves) - 1, 1)), lw=0.8)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xscale("log")
    ax.set_title("(a) posterior-mode curve at each refit (dark: early, light: late)")
    ax.set_xlabel("lag in days (log scale)")
    ax.set_ylabel("coefficient (scaled return per scaled return)")
    ax = axes[1]
    for r in results[:3]:
        pnl = r["pnl"] - 0.0005 * r["traded"]
        ax.plot(close.index[r["time"]], np.cumsum(pnl) / pnl.std() / np.sqrt(252), label=f"badlm {r['rule']}")
    for r in results[3:]:
        if r["rule"].startswith("gp"):
            pnl = r["pnl"] - 0.0005 * r["traded"]
            ax.plot(close.index[r["time"]], np.cumsum(pnl) / pnl.std() / np.sqrt(252), ls="--", lw=1, label=r["signal"])
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("(b) cumulative return net of 5 bp, unit annual volatility")
    ax = axes[2]
    ax.plot(lags, np.cumsum(curves[-1]))
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xscale("log")
    ax.set_title("(c) cumulative curve: expected scaled return from a unit shock j days ago")
    ax.set_xlabel("lag in days (log scale)")
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, f"trading_{tag}.png"), dpi=130)

    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"trading_{tag}.json"), "w") as f:
        json.dump({"args": vars(args), "n_assets": int(df["asset"].nunique()), "oos_r2": oos_r2,
                   "ic": ic(rows["pred"], rows["y"]), "n_refits": len(wf.refit_times),
                   "results": strip(results), "by_kind": by_kind, "deflated_sharpe": dsr,
                   "yearly_sharpe": {int(k): float(v) for k, v in ys.items()},
                   "final_curve": curves[-1].tolist()}, f, indent=1)
    print("written", flush=True)
