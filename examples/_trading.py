"""Shared pieces of the trading experiments: stacked rows of volatility-scaled
returns, benchmark signals, cost-aware evaluation and result tables."""
from __future__ import annotations

import numpy as np
import pandas as pd

from _markets import ewma_vol, log_returns

from badlm.backtest import performance, trade

COSTS = (0.0, 0.0001, 0.0002, 0.0005, 0.001)


def stacked_rows(close, vol_halflife=30, clip=5.0, kinds=None, use=("equity", "fx", "commodity"),
                 volume=None, min_days=600):
    """Rows (asset, t): ``x`` the scaled return on day t, ``y`` the scaled return
    on day t + 1, the raw next return and the volatility that scales it, the
    panel position of the target day and, if ``volume`` is given, the log of
    volume relative to its trailing 60-day mean on day t."""
    r = log_returns(close)
    vol = ewma_vol(r, halflife=vol_halflife)
    rs = (r / vol).clip(-clip, clip)
    parts = []
    for a in close:
        if kinds is not None and kinds[a] not in use:
            continue
        s = rs[a].dropna()
        if s.size < min_days:
            continue
        pos = close.index.get_indexer(s.index)
        part = {"asset": a, "x": s.to_numpy()[:-1], "y": s.to_numpy()[1:],
                "raw": r[a].reindex(s.index).to_numpy()[1:],
                "vol": vol[a].reindex(s.index).to_numpy()[1:], "time": pos[1:]}
        if volume is not None and a in volume:
            v = volume[a].reindex(s.index).astype(float)
            rel = np.log(v.clip(lower=1.0)) - np.log(v.clip(lower=1.0)).rolling(60, min_periods=20).mean()
            part["relvol"] = rel.fillna(0.0).to_numpy()[:-1]
        parts.append(pd.DataFrame(part))
    return pd.concat(parts, ignore_index=True)


def benchmark_features(df, windows=(5, 21, 63, 252)):
    """Trailing sums of the scaled exposure as of the row's x date, per asset."""
    out = {}
    for w in windows:
        out[f"sum{w}"] = df.groupby("asset")["x"].transform(lambda s: s.rolling(w).sum())
    return pd.DataFrame(out).fillna(0.0)


def evaluate(name, aim_scaled, rows, costs=COSTS, rules=(("aim", {}), ("gp", {"speed": 0.25}), ("gp", {"speed": 0.1}))):
    """Positions from an aim in scaled units (converted to notional by 1/vol)
    under several trading rules; returns one result dict per rule."""
    res = []
    var = np.ones_like(aim_scaled)
    for rule, kw in rules:
        tr = trade(aim_scaled, var, rows["asset"].to_numpy(), rows["time"].to_numpy(),
                   vol=rows["vol"].to_numpy(), rule=rule, **kw)
        perf = performance(tr["position"], tr["trade"], rows["raw"].to_numpy(), rows["time"].to_numpy(), costs=costs)
        label = rule + ("" if not kw else " " + " ".join(f"{v}" for v in kw.values()))
        res.append({"signal": name, "rule": label, "sharpe_gross": perf["sharpe_gross"],
                    "turnover": perf["turnover_annual"], "breakeven_bp": 1e4 * perf["breakeven_cost"],
                    **{f"net_{int(round(c * 1e4))}bp": v for c, v in perf["sharpe_net"].items()},
                    "pnl": perf["pnl_gross"], "traded": perf["traded"], "time": perf["time"]})
    return res


def ic(pred, y):
    return float(np.corrcoef(np.asarray(pred, float), np.asarray(y, float))[0, 1])


def print_table(results, title):
    print(f"\n{title}")
    print("| Signal | Rule | Gross SR | Turnover/yr | Break-even (bp) | Net SR 1bp | 2bp | 5bp | 10bp |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in results:
        print(f"| {r['signal']} | {r['rule']} | {r['sharpe_gross']:.2f} | {r['turnover']:.0f} | {r['breakeven_bp']:.1f} | "
              f"{r['net_1bp']:.2f} | {r['net_2bp']:.2f} | {r['net_5bp']:.2f} | {r['net_10bp']:.2f} |")


def strip(results):
    return [{k: v for k, v in r.items() if k not in ("pnl", "traded", "time")} for r in results]


def yearly_sharpe(result, index, cost=0.0005):
    pnl = result["pnl"] - cost * result["traded"]
    t = index[result["time"]]
    g = pd.Series(pnl, index=t).groupby(t.year).agg(["mean", "std"])
    return (g["mean"] / g["std"] * np.sqrt(252)).round(2)
