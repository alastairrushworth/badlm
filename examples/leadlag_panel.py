"""Pillar A: a rank-R lag-curve matrix across assets.

Every asset's next-day scaled return responds to a few driver series through
shared adaptive lag curves and asset-specific loadings shrunk towards one
(``badlm.panel``). Drivers: the asset's own past, and the past of the
equal-risk baskets of ``factor_curves.py`` (US market, UK market, dollar,
commodities). The own-past curve is a surface over the asset's own
volatility state; the basket drivers have plain curves. Three fits are compared walk-forward on the same rows:

(a) own past only, loadings fixed at one: the pooled single-series surface;
(b) own past only with hierarchical loadings;
(c) own past plus the basket drivers with hierarchical loadings: the
    lead-lag matrix of rank five.

Run: python examples/leadlag_panel.py [--final] [--universe equity|fx|all]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _markets import ewma_vol, log_returns  # noqa: E402
from _panel import HOLDOUT_START, load_wide_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows, strip  # noqa: E402
from factor_curves import DROP, basket_returns, baskets  # noqa: E402

from badlm.backtest import subperiod_sharpe  # noqa: E402
from badlm.panel import build_panel_design, panel_walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")


def report(name, wf, df, close, breaks):
    rows = df.iloc[wf["row"]].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf["mean"], wf["var"]
    res = evaluate(name, (rows["pred"] / rows["pvar"]).to_numpy(), rows, rules=(("aim", {}), ("gp", {"speed": 0.05})))
    slow = res[1]
    sub = subperiod_sharpe(slow["pnl"] - 0.0005 * slow["traded"], slow["time"], breaks)
    print(f"{name}: IC {ic(rows['pred'], rows['y']):.4f}; sub-period net Sharpe "
          f"{[(close.index[int(lo)].year, round(s, 2)) for lo, hi, s, n in sub]}", flush=True)
    return res, rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--universe", default="all")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--loading-scale", type=float, default=0.5)
    args = ap.parse_args()
    panel = load_wide_panel()
    close, kind = panel["close"], panel["kind"]
    close = close.drop(columns=[c for c in DROP if c in close])
    if not args.final:
        close = close.loc[:HOLDOUT_START - pd.Timedelta(days=1)]
    memb = baskets(close, kind)
    drivers_b = ["us_market", "uk_market", "dollar", "commodity"]
    bret = pd.DataFrame({k: basket_returns(close, memb[k]) for k in drivers_b})
    bvol = ewma_vol(bret)
    bscaled = (bret / bvol).clip(-5, 5).fillna(0.0)
    use = {"equity": ("equity_us", "equity_uk", "index"), "fx": ("fx",),
           "all": ("equity_us", "equity_uk", "index", "fx", "commodity")}[args.universe]
    df = stacked_rows(close, kinds=kind, use=use)
    med = df.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    df["volstate"] = np.log(df["vol"] / med).fillna(0.0).clip(-2, 2)
    xdate = close.index[df["time"] - 1]
    for b in drivers_b:
        df[f"x_{b}"] = bscaled[b].reindex(xdate).to_numpy()
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    breaks = [int(close.index.get_indexer([pd.Timestamp(y)], method="bfill")[0]) for y in ("2004-01-01", "2012-01-01")]
    pri = Priors(tau_level_scale=0.001, temper=10)
    print(f"{df['asset'].nunique()} assets, {len(df)} rows, drivers: own + {drivers_b}", flush=True)
    t_all = df["time"].to_numpy(); a_all = df["asset"].to_numpy(); y_all = df["y"].to_numpy()
    results, out = [], {}

    t0 = time.time()
    d1 = build_panel_design([df["x"].to_numpy()], y_all, t_all, a_all, nlag=250, k=args.k,
                            modifier=df["volstate"].to_numpy(), modifier_k=5)
    wf = panel_walk_forward(d1, start=start, step=args.step, priors=pri, pool=True)
    res, rows = report("(a) own past, pooled", wf, df, close, breaks); results += res
    out["a"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    wf = panel_walk_forward(d1, start=start, step=args.step, priors=pri, loading_scale=args.loading_scale)
    res, rows = report("(b) own past, hierarchical loadings", wf, df, close, breaks); results += res
    ld = wf["refits"][-1]["loading"][:, 0]
    print("   final loadings on own past: mean %.2f sd %.2f, range %.2f to %.2f" % (ld.mean(), ld.std(), ld.min(), ld.max()))
    out["b"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"]),
                "loadings": dict(zip(d1.names, ld.round(3).tolist()))}
    drivers = [df["x"].to_numpy()] + [df[f"x_{b}"].to_numpy() for b in drivers_b]
    d2 = build_panel_design(drivers, y_all, t_all, a_all, nlag=250, k=args.k,
                            modifier=df["volstate"].to_numpy(), modifier_k=5, modifier_drivers=[0])
    wf = panel_walk_forward(d2, start=start, step=args.step, priors=pri, loading_scale=args.loading_scale)
    res, rows = report("(c) own + basket drivers, hierarchical loadings", wf, df, close, breaks); results += res
    last = wf["refits"][-1]
    print("   loading sd per driver:", dict(zip(["own"] + drivers_b, np.round(last["loading_sd"], 3))))
    for r, name in enumerate(["own"] + drivers_b):
        c = last["curves"][r][2] if r == 0 else last["curves"][r]     # own: surface at the middle state
        print(f"   {name:10s} curve (mid state): lag1 {c[0]:+.4f} lag5 {c[4]:+.4f} lag21 {c[20]:+.4f} sum {c.sum():+.3f}; "
              f"loading mean {last['loading'][:, r].mean():.2f} sd {last['loading'][:, r].std():.2f}")
    out["c"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"]),
                "loading_sd": dict(zip(["own"] + drivers_b, np.round(last["loading_sd"], 4).tolist())),
                "loadings": {name: dict(zip(d2.names, last["loading"][:, r].round(3).tolist())) for r, name in enumerate(["own"] + drivers_b)}}
    print(f"three panel fits in {time.time() - t0:.0f}s")
    print_table(results, f"{args.universe}: lead-lag panel, out of sample (development period)")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"leadlag_{args.universe}{'_final' if args.final else ''}.json"), "w") as f:
        json.dump({"args": vars(args), **out}, f, indent=1)
    print("written")
