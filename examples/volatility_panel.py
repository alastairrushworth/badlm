"""Pillar C: a lag-curve model of volatility across assets, and sizing with it.

The target is the log absolute return of each asset relative to its slow
volatility level (a 250-day EWMA), so that assets with different volatility
share a model; the drivers are the asset's own past of that series and the
cross-sectional mean of it over the asset's class (the volatility spillover
of Diebold and Yilmaz, with the lag structure estimated). Three forecasters
are compared on the QLIKE loss of the implied variance against the next
day's squared return: the EWMA used for sizing so far, a HAR-type rule on
the same series (means over 1, 5 and 22 days, pooled least squares), the
pooled lag curve, the hierarchical panel and the panel with the spillover
driver. The best forecast then replaces the EWMA in sizing the asset-level
volatility-surface strategy.

Run: python examples/volatility_panel.py [--final]
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
from _trading import evaluate, ic, print_table, stacked_rows  # noqa: E402
from factor_curves import DROP  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.backtest import walk_forward  # noqa: E402
from badlm.panel import build_panel_design, panel_walk_forward  # noqa: E402
from badlm.priors import Priors  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
resdir = os.path.join(here, "..", "paper", "results")
CLASS = {"equity_us": "us", "equity_uk": "uk", "index": "us", "fx": "fx", "commodity": "cmd", "factor_etf": "us"}


def qlike(r2, var):
    """Patton's QLIKE loss; squared returns are floored at a hundredth of
    the forecast variance so that zero-return days (stale prices) do not
    dominate."""
    z = np.maximum(r2, 0.01 * var) / var
    return float(np.mean(z - np.log(z) - 1.0))


def vol_rows(close, kind):
    """Rows (asset, t) with the relative log absolute return on day t (x) and
    day t + 1 (y), the slow volatility level, the class mean of x, and the
    squared raw return of day t + 1."""
    r = log_returns(close)
    slow = r.ewm(halflife=250, min_periods=120).std().shift(1)
    a = np.log(r.abs() + 0.05 * slow) - np.log(slow)          # relative log absolute return
    cls = pd.Series({c: CLASS[kind[c]] for c in close})
    mkt = pd.DataFrame({g: a[cls.index[cls == g]].mean(axis=1) for g in cls.unique()})
    parts = []
    for c in close:
        s = a[c].dropna()
        if s.size < 600:
            continue
        pos = close.index.get_indexer(s.index)
        parts.append(pd.DataFrame({"asset": c, "x": s.to_numpy()[:-1], "y": s.to_numpy()[1:],
                                   "x_mkt": mkt[cls[c]].reindex(s.index).fillna(0.0).to_numpy()[:-1],
                                   "slow": slow[c].reindex(s.index).to_numpy()[1:],
                                   "r2": (r[c].reindex(s.index).to_numpy()[1:]) ** 2,
                                   "ewma30": ewma_vol(r)[c].reindex(s.index).to_numpy()[1:],
                                   "time": pos[1:]}))
    return pd.concat(parts, ignore_index=True)


def har_forecast(df, start, step):
    """Pooled least-squares HAR on the relative log absolute return, refitted
    at the same cadence as the lag-curve models."""
    g = df.groupby("asset")["x"]
    F = np.c_[np.ones(len(df)), df["x"], g.transform(lambda s: s.rolling(5).mean()), g.transform(lambda s: s.rolling(22).mean())]
    ok = np.isfinite(F).all(axis=1)
    pred = np.full(len(df), np.nan)
    t = df["time"].to_numpy()
    T = start
    while T < t.max():
        train = ok & (t <= T)
        test = ok & (t > T) & (t <= T + step)
        if train.sum() > 100 and test.any():
            beta, *_ = np.linalg.lstsq(F[train], df["y"].to_numpy()[train], rcond=None)
            pred[test] = F[test] @ beta
        T += step
    return pred


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--step", type=int, default=63)
    args = ap.parse_args()
    panel = load_wide_panel()
    close, kind = panel["close"], panel["kind"]
    close = close.drop(columns=[c for c in DROP if c in close] + [c for c in close if kind[c] == "factor_etf"])
    if not args.final:
        close = close.loc[:HOLDOUT_START - pd.Timedelta(days=1)]
    df = vol_rows(close, kind)
    start = int(close.index.get_indexer([pd.Timestamp(f"{args.start}-01-01")], method="bfill")[0])
    pri = Priors(tau_level_scale=0.01, temper=10)
    print(f"{df['asset'].nunique()} assets, {len(df)} rows", flush=True)
    t_all, a_all, y_all = df["time"].to_numpy(), df["asset"].to_numpy(), df["y"].to_numpy()
    forecasts = {}

    # a variance forecast from a log|r| prediction m with residual variance v: E|r| = slow exp(m + v/2) * c0,
    # calibrated by the ratio of mean r^2 to mean (slow exp(m))^2 on the training rows (done implicitly by
    # a pooled scale factor fitted on the first ten years)
    def variance_from(pred, var, rows):
        base = (rows["slow"] * np.exp(pred + 0.5 * var)) ** 2
        early = rows["time"] < start + 252 * 10
        scale = np.mean(rows.loc[early, "r2"]) / np.mean(base[early])
        return base * scale

    t0 = time.time()
    d1 = build_panel_design([df["x"].to_numpy()], y_all, t_all, a_all, nlag=250, k=25)
    for name, kw in (("pooled curve", {"pool": True}), ("hierarchical panel", {})):
        wf = panel_walk_forward(d1, start=start, step=args.step, priors=pri, **kw)
        rows = df.iloc[wf["row"]].reset_index(drop=True)
        forecasts[name] = (rows, wf["mean"], wf["var"])
        print(f"{name}: IC on log|r| {ic(wf['mean'], wf['y']):.3f} ({time.time() - t0:.0f}s)", flush=True)
    d2 = build_panel_design([df["x"].to_numpy(), df["x_mkt"].to_numpy()], y_all, t_all, a_all, nlag=250, k=25)
    wf = panel_walk_forward(d2, start=start, step=args.step, priors=pri)
    rows = df.iloc[wf["row"]].reset_index(drop=True)
    forecasts["panel + class spillover"] = (rows, wf["mean"], wf["var"])
    last = wf["refits"][-1]
    print(f"panel + class spillover: IC on log|r| {ic(wf['mean'], wf['y']):.3f}; own curve lag1 {last['curves'][0][0]:+.3f} "
          f"lag5 {last['curves'][0][4]:+.3f} lag22 {last['curves'][0][21]:+.3f} sum {last['curves'][0].sum():+.2f}; "
          f"spillover curve lag1 {last['curves'][1][0]:+.3f} lag5 {last['curves'][1][4]:+.3f} sum {last['curves'][1].sum():+.2f}", flush=True)
    har = har_forecast(df, start, args.step)

    # QLIKE on common rows
    base_rows = forecasts["pooled curve"][0]
    key = base_rows[["asset", "time"]]
    print("\n| Forecast | QLIKE | relative to EWMA(30) |")
    print("|---|---|---|")
    out = {}
    ew = qlike(base_rows["r2"].to_numpy(), base_rows["ewma30"].to_numpy() ** 2)
    out["ewma30"] = ew
    print(f"| EWMA, 30-day half-life | {ew:.4f} | 1.000 |")
    hv = pd.Series(har, index=pd.MultiIndex.from_frame(df[["asset", "time"]])).reindex(pd.MultiIndex.from_frame(key)).to_numpy()
    okh = np.isfinite(hv)
    hvar = variance_from(hv[okh], np.var(base_rows["y"].to_numpy()[okh] - hv[okh]), base_rows[okh].reset_index(drop=True))
    q = qlike(base_rows["r2"].to_numpy()[okh], hvar.to_numpy())
    out["har"] = q
    print(f"| HAR, pooled least squares | {q:.4f} | {q / ew:.3f} |")
    for name, (rows, mean, var) in forecasts.items():
        v = variance_from(mean, var, rows)
        q = qlike(rows["r2"].to_numpy(), v.to_numpy())
        out[name] = q
        print(f"| {name} | {q:.4f} | {q / ew:.3f} |")
    best_name = min(forecasts, key=lambda n: out[n])
    print(f"best lag-curve forecaster: {best_name}")

    # sizing the asset-level volatility-surface strategy with the best forecast
    rows_v, mean_v, var_v = forecasts[best_name]
    vol_fc = pd.Series(np.sqrt(variance_from(mean_v, var_v, rows_v).to_numpy()),
                       index=pd.MultiIndex.from_frame(rows_v[["asset", "time"]]))
    dfs = stacked_rows(close, kinds=kind, use=("equity_us", "equity_uk", "index", "fx", "commodity"))
    med = dfs.groupby("asset")["vol"].transform(lambda s: s.expanding(min_periods=60).median())
    dfs["volstate"] = np.log(dfs["vol"] / med).fillna(0.0).clip(-2, 2)
    wfs = walk_forward(BADLM(nlag=250, k=25, modifier_k=5, priors=Priors(tau_level_scale=0.001, temper=10)),
                       dfs["x"].to_numpy(), dfs["y"].to_numpy(), dfs["time"].to_numpy(), groups=dfs["asset"].to_numpy(),
                       start=start, step=args.step, modifier=dfs["volstate"].to_numpy())
    rs = dfs.iloc[wfs.row].reset_index(drop=True)
    rs["pred"], rs["pvar"] = wfs.mean, wfs.var
    fc = vol_fc.reindex(pd.MultiIndex.from_frame(rs[["asset", "time"]])).to_numpy()
    ok = np.isfinite(fc)
    rs = rs[ok].reset_index(drop=True); fc = fc[ok]
    results = evaluate("surface, EWMA sizing", (rs["pred"] / rs["pvar"]).to_numpy(), rs, rules=(("gp", {"speed": 0.05}),))
    rs2 = rs.copy(); rs2["vol"] = fc
    results += evaluate(f"surface, sized with {best_name}", (rs2["pred"] / rs2["pvar"]).to_numpy(), rs2, rules=(("gp", {"speed": 0.05}),))
    rs3 = rs.copy(); rs3["vol"] = np.sqrt(rs["vol"] * fc)
    results += evaluate("surface, geometric mean of the two", (rs3["pred"] / rs3["pvar"]).to_numpy(), rs3, rules=(("gp", {"speed": 0.05}),))
    print_table(results, "Sizing the asset-level surface strategy (development period)")
    os.makedirs(resdir, exist_ok=True)
    with open(os.path.join(resdir, f"volatility_panel{'_final' if args.final else ''}.json"), "w") as f:
        json.dump({"args": vars(args), "qlike": out, "best": best_name,
                   "sizing": [{k: v for k, v in r.items() if k not in ("pnl", "traded", "time")} for r in results]}, f, indent=1)
    print("written")
