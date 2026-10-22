"""Experiment 2: a lag curve that varies with the volatility state.

The curve becomes a surface over lag and a modifier ``z`` with the adaptive
penalty along lag and a penalised-complexity prior on roughness along ``z``
(``BADLM(modifier_k=...)``). The modifier is the asset's own volatility
state, the log of its ex-ante volatility relative to its expanding median,
after Nagel (2012) and Daniel and Moskowitz (2016). The static pooled curve
is fitted on the same rows for comparison.

Run: python examples/trading_modulated.py [--universe fx] [options]
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


def run_variant(name, model, df, start, step, modifier=None):
    t0 = time.time()
    wf = walk_forward(model, df["x"].to_numpy(), df["y"].to_numpy(), df["time"].to_numpy(),
                      groups=df["asset"].to_numpy(), start=start, step=step, modifier=modifier)
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
    ap.add_argument("--modifier-k", type=int, default=5)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--temper", type=float, default=1.0)
    args = ap.parse_args()
    panel = load_panel()
    close, kinds = panel["close"], panel["kind"]
    use = ("equity", "fx", "commodity") if args.universe == "all" else (args.universe,)
    df = stacked_rows(close, kinds=kinds, use=use)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    pri = Priors(tau_level_scale=args.scale, temper=args.temper)
    # volatility state known on the x date: the vol that scales y is the ex-ante vol of day t + 1,
    # formed from returns up to t
    med = df.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    df["volstate"] = np.log(df["vol"] / med).fillna(0.0).clip(-2, 2)
    print(f"{df['asset'].nunique()} assets, {len(df)} rows", flush=True)

    out, results = {}, []
    res, wf, rows = run_variant("static curve", BADLM(nlag=args.nlag, k=args.k, priors=pri), df, start, args.step)
    out["static"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"]), "curve": wf.curves[-1][0].tolist()}
    results += res
    res, wfv, rows = run_variant("volatility-modulated", BADLM(nlag=args.nlag, k=args.k, modifier_k=args.modifier_k, priors=pri),
                                 df, start, args.step, modifier=df["volstate"].to_numpy())
    surf = wfv.curves[-1][0]
    print("final surface at modifier grid (low to high volatility), sum over lags:",
          " ".join(f"{s.sum():+.3f}" for s in surf), "| lag 1-5 mean:", " ".join(f"{s[:5].mean():+.4f}" for s in surf))
    out["volatility"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"]), "surface": surf.tolist()}
    results += res
    print_table(results, f"{args.universe}: modulated curves, out of sample")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"trading_modulated_{args.universe}.json"), "w") as f:
        json.dump({"args": vars(args), **out}, f, indent=1)
    print("written")
