"""Pillar B: lag curves on factor baskets, with state surfaces.

Assets are grouped into equal-risk baskets (US market, UK market, US sectors,
the dollar against all currencies, commodities, the index ETFs); each basket's
daily return is the mean of its members' volatility-scaled returns. A pooled
adaptive curve of the next day's basket return on its own past, with the
volatility-state surface, is estimated walk-forward and compared with the
asset-level surface on the same members. A second variant adds the basket's
average pairwise correlation as a second, additive state surface (one
modifier per exposure), and a third adds the market basket's past as a
second exposure.

Development runs stop at ``HOLDOUT_START``; ``--final`` includes the held-out
years. Run: python examples/factor_curves.py [--final]
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
from _panel import HOLDOUT_START, SECTOR, load_wide_panel  # noqa: E402
from _trading import evaluate, ic, print_table, stacked_rows, strip  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import subperiod_sharpe, walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")
DROP = ("RRC", "JD.L", "SHLD")


def baskets(close, kind):
    """Membership of the baskets."""
    b = {"us_market": [c for c in close if kind[c] == "equity_us"],
         "uk_market": [c for c in close if kind[c] == "equity_uk"],
         "dollar": [c for c in close if kind[c] == "fx"],
         "commodity": [c for c in close if kind[c] == "commodity"],
         "index_etf": [c for c in close if kind[c] == "index"]}
    for sec in ("tech", "fin", "energy", "health", "cons"):
        b[f"us_{sec}"] = [c for c in close if kind[c] == "equity_us" and SECTOR.get(c) == sec]
    return {k: v for k, v in b.items() if len(v) >= 3}


def basket_returns(close, members, halflife=30):
    """Mean of the members' scaled returns, as a daily return series of the
    basket (so that a basket is traded like an asset with unit-vol members)."""
    r = log_returns(close[members])
    vol = ewma_vol(r, halflife=halflife)
    rs = (r / vol).clip(-5, 5)
    n = rs.notna().sum(axis=1)
    return rs.mean(axis=1).where(n >= 3) / 100.0        # in 'per cent of a unit-vol position'


def state_vars(df, close_basket_members, rs_all):
    """Own volatility state (as before) and the basket's average pairwise
    correlation over the past 60 days, known on the x date."""
    med = df.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    df["volstate"] = np.log(df["vol"] / med).fillna(0.0).clip(-2, 2)
    corr = {}
    for name, members in close_basket_members.items():
        sub = rs_all[members]
        # mean pairwise correlation from the variance of the row mean: var(mean) = v (1 + (n - 1) rho) / n
        v = sub.rolling(60, min_periods=40).var()
        rowmean_var = sub.mean(axis=1).rolling(60, min_periods=40).var()
        n = sub.notna().sum(axis=1).clip(lower=1)
        corr[name] = ((n * rowmean_var - v.mean(axis=1) / 1.0) / (n - 1).clip(lower=1) / v.mean(axis=1)).clip(-1, 1)
    corr = pd.DataFrame(corr)
    df["corrstate"] = [corr.at[t, a] if a in corr else np.nan for a, t in zip(df["asset"], df["time_dt"])]
    df["corrstate"] = df["corrstate"].fillna(0.0).clip(-0.2, 1.0)
    return df


def run(name, model, df, X, start, step, modifier=None):
    t0 = time.time()
    wf = walk_forward(model, X, df["y"].to_numpy(), df["time"].to_numpy(), groups=df["asset"].to_numpy(),
                      start=start, step=step, modifier=modifier)
    rows = df.iloc[wf.row].reset_index(drop=True)
    rows["pred"], rows["pvar"] = wf.mean, wf.var
    res = evaluate(name, (rows["pred"] / rows["pvar"]).to_numpy(), rows, rules=(("aim", {}), ("gp", {"speed": 0.05})))
    print(f"{name}: {len(wf.refit_times)} refits in {time.time() - t0:.0f}s, IC {ic(rows['pred'], rows['y']):.4f}", flush=True)
    return res, rows, wf


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    ap.add_argument("--k", type=int, default=25)
    args = ap.parse_args()
    panel = load_wide_panel()
    close, kind = panel["close"], panel["kind"]
    close = close.drop(columns=[c for c in DROP if c in close])
    if not args.final:
        close = close.loc[:HOLDOUT_START - pd.Timedelta(days=1)]
    memb = baskets(close, kind)
    print("baskets:", {k: len(v) for k, v in memb.items()}, flush=True)
    r_all = log_returns(close)
    rs_all = (r_all / ewma_vol(r_all)).clip(-5, 5)
    bret = pd.DataFrame({k: basket_returns(close, v) for k, v in memb.items()})
    bclose = np.exp(bret.fillna(0.0).cumsum()).where(bret.notna().cummax())
    bkind = {k: "basket" for k in memb}
    dfb = stacked_rows(bclose, kinds=bkind, use=("basket",))
    dfb["time_dt"] = close.index[dfb["time"]]
    dfb = state_vars(dfb, memb, rs_all)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    pri = Priors(tau_level_scale=0.001, temper=10)
    breaks = [int(close.index.get_indexer([pd.Timestamp(y)], method="bfill")[0]) for y in ("2004-01-01", "2012-01-01")]
    out, results = {}, []

    xb = dfb["x"].to_numpy()
    res, rows, wf = run("baskets: static curve", BADLM(nlag=250, k=args.k, priors=pri), dfb, xb, start, args.step)
    results += res; out["basket_static"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    res, rows, wf = run("baskets: volatility surface", BADLM(nlag=250, k=args.k, modifier_k=5, priors=pri), dfb, xb,
                        start, args.step, modifier=dfb["volstate"].to_numpy())
    results += res; out["basket_vol"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    best = res[1]
    print("  sub-period Sharpe (slow, 5 bp):", [(close.index[int(lo)].year, round(s, 2)) for lo, hi, s, n in
          subperiod_sharpe(best["pnl"] - 0.0005 * best["traded"], best["time"], breaks)])
    # additive surfaces: own exposure twice, one with the volatility state and one with the correlation state
    X2 = np.c_[xb, xb]
    Z2 = np.c_[dfb["volstate"].to_numpy(), dfb["corrstate"].to_numpy()]
    res, rows, wf = run("baskets: volatility + correlation surfaces", BADLM(nlag=[250, 250], k=[args.k, args.k], modifier_k=5, priors=pri),
                        dfb, X2, start, args.step, modifier=Z2)
    results += res; out["basket_vol_corr"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    # cross-basket: the US market basket's past as a second exposure for every basket
    mkt = bclose["us_market"]
    mkt_rs = (log_returns(mkt) / ewma_vol(log_returns(mkt))).clip(-5, 5)
    dfb["x_mkt"] = mkt_rs.reindex(close.index[dfb["time"] - 1]).to_numpy()
    dfb["x_mkt"] = dfb["x_mkt"].fillna(0.0)
    Xm = np.c_[xb, dfb["x_mkt"].to_numpy()]
    res, rows, wf = run("baskets: own + US market past, volatility surface",
                        BADLM(nlag=[250, 250], k=[args.k, args.k], modifier_k=5, priors=pri), dfb, Xm, start, args.step,
                        modifier=dfb["volstate"].to_numpy())
    results += res; out["basket_vol_mkt"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    cm = wf.curves[-1][1]
    print("  US-market-past curve at the average state: lag1 %.4f lag5 %.4f lag21 %.4f sum %.3f" % (cm[2][0], cm[2][4], cm[2][20], cm[2].sum()))
    print_table(results, "Baskets, out of sample (development period)")

    # asset-level surface on the members of the same baskets, for comparison
    members = sorted({m for v in memb.values() for m in v})
    dfa = stacked_rows(close[members], kinds={m: "x" for m in members}, use=("x",))
    med = dfa.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    dfa["volstate"] = np.log(dfa["vol"] / med).fillna(0.0).clip(-2, 2)
    res, rows, wf = run("assets: volatility surface", BADLM(nlag=250, k=args.k, modifier_k=5, priors=pri), dfa,
                        dfa["x"].to_numpy(), start, args.step, modifier=dfa["volstate"].to_numpy())
    out["asset_vol"] = {"results": strip(res), "ic": ic(rows["pred"], rows["y"])}
    best = res[1]
    print("  sub-period Sharpe (slow, 5 bp):", [(close.index[int(lo)].year, round(s, 2)) for lo, hi, s, n in
          subperiod_sharpe(best["pnl"] - 0.0005 * best["traded"], best["time"], breaks)])
    print_table(res, "Members traded individually, out of sample (development period)")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"factor_curves{'_final' if args.final else ''}.json"), "w") as f:
        json.dump({"args": vars(args), "baskets": {k: v for k, v in memb.items()}, **out}, f, indent=1)
    print("written")
