"""Experiment 2: several exposures with their own adaptive lag curves.

(a) Asymmetry: positive and negative past scaled returns as two exposures, so
    the response to gains and to losses can have different lag curves.
(b) Volume interaction (equities): the scaled return and the scaled return
    multiplied by relative volume (log volume over its trailing 60-day mean)
    as two exposures, after Llorente, Michaely, Saar and Wang (2002): moves on
    high volume tend to continue, moves on low volume to revert.

Each variant is evaluated walk-forward against the single-exposure curve on
identical rows. Run: python examples/trading_interactions.py [--universe fx|equity|all] [options]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _markets import load_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows, strip  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")


def run_variant(name, model, X, df, start, step):
    t0 = time.time()
    wf = walk_forward(model, X, df["y"].to_numpy(), df["time"].to_numpy(), groups=df["asset"].to_numpy(),
                      start=start, step=step)
    rows = df.iloc[wf.row].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf.mean, wf.var
    res = evaluate(name, (rows["pred"] / rows["pvar"]).to_numpy(), rows, rules=(("aim", {}), ("gp", {"speed": 0.1})))
    print(f"{name}: {len(wf.refit_times)} refits in {time.time() - t0:.0f}s, IC {ic(rows['pred'], rows['y']):.4f}", flush=True)
    return res, wf, rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="fx")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--nlag", type=int, default=250)
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--temper", type=float, default=1.0)
    args = ap.parse_args()
    panel = load_panel()
    close, kinds = panel["close"], panel["kind"]
    use = ("equity", "fx", "commodity") if args.universe == "all" else (args.universe,)
    df = stacked_rows(close, kinds=kinds, use=use, volume=panel["volume"])
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    pri = Priors(tau_level_scale=args.scale, temper=args.temper)
    print(f"{df['asset'].nunique()} assets, {len(df)} rows", flush=True)

    out = {}
    x = df["x"].to_numpy()
    res1, wf1, rows1 = run_variant("single curve", BADLM(nlag=args.nlag, k=args.k, priors=pri), x, df, start, args.step)
    out["single"] = {"results": strip(res1), "ic": ic(rows1["pred"], rows1["y"]), "curve": wf1.curves[-1][0].tolist()}
    Xa = np.c_[np.maximum(x, 0.0), np.minimum(x, 0.0)]
    res2, wf2, rows2 = run_variant("gains + losses", BADLM(nlag=[args.nlag] * 2, k=[args.k] * 2, priors=pri),
                                   Xa, df, start, args.step)
    out["asymmetric"] = {"results": strip(res2), "ic": ic(rows2["pred"], rows2["y"]),
                         "curve_gains": wf2.curves[-1][0].tolist(), "curve_losses": wf2.curves[-1][1].tolist()}
    cg, cl = wf2.curves[-1]
    print("final curves, gains vs losses at lags 1, 5, 21, 63:",
          " ".join(f"{cg[i]:+.4f}/{cl[i]:+.4f}" for i in (0, 4, 20, 62)), f"sums {cg.sum():+.3f}/{cl.sum():+.3f}")
    results = res1 + res2
    if "relvol" in df:
        has = df["relvol"].notna().to_numpy()
        dv = df[has].reset_index(drop=True)
        xv = dv["x"].to_numpy()
        Xv = np.c_[xv, xv * dv["relvol"].to_numpy()]
        res3, wf3, rows3 = run_variant("return + return x volume",
                                       BADLM(nlag=[args.nlag] * 2, k=[args.k] * 2, priors=pri), Xv, dv, start, args.step)
        res0, wf0, rows0 = run_variant("single curve (same rows)", BADLM(nlag=args.nlag, k=args.k, priors=pri),
                                       xv, dv, start, args.step)
        c0, cv = wf3.curves[-1]
        print("final curves, return vs return x volume at lags 1, 5, 21, 63:",
              " ".join(f"{c0[i]:+.4f}/{cv[i]:+.4f}" for i in (0, 4, 20, 62)))
        out["volume"] = {"results": strip(res3 + res0), "ic": ic(rows3["pred"], rows3["y"]),
                         "ic_single": ic(rows0["pred"], rows0["y"]),
                         "curve_return": c0.tolist(), "curve_volume": cv.tolist()}
        results += res3 + res0
    print_table(results, f"{args.universe}: multi-exposure curves, out of sample")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"trading_interactions_{args.universe}.json"), "w") as f:
        json.dump({"args": vars(args), **out}, f, indent=1)
    print("written")
