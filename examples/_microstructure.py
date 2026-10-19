"""Turn Lean tick files (trades and quotes) into per-second microstructure series.

The sample data shipped with the QuantConnect Lean repository hold every trade
and quote update for SPY, BAC, AIG and IBM on 7 to 11 October 2013. For each
regular-session second this module computes the national best bid and offer
(the best bid and ask over the latest quote of each exchange), the midquote
return in basis points, the number and volume of buyer- and seller-initiated
trades (trades above the prevailing midquote are buys, below are sells, and
trades at the midquote are signed by the tick rule), and the quoted spread.

Requires pandas. Results are cached in ``examples/data/micro_<symbol>.npz``.
"""
from __future__ import annotations

import os
import urllib.request

import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
DAYS = ("20131007", "20131008", "20131009", "20131010", "20131011")
OPEN_MS, CLOSE_MS = 34_200_000, 57_600_000           # 09:30 to 16:00
N_SEC = (CLOSE_MS - OPEN_MS) // 1000
URL = "https://raw.githubusercontent.com/QuantConnect/Lean/master/Data/equity/usa/tick/{sym}/{day}_{kind}.zip"


def _fetch(sym, day, kind):
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, f"{sym}_{day}_{kind}.zip")
    if not os.path.exists(path):
        print(f"downloading {sym} {day} {kind}")
        urllib.request.urlretrieve(URL.format(sym=sym, day=day, kind=kind), path)
    return path


def process_day(sym, day):
    import pandas as pd

    trades = pd.read_csv(_fetch(sym, day, "trade"), header=None,
                         names=["ms", "price", "size", "exch", "cond", "susp"])
    quotes = pd.read_csv(_fetch(sym, day, "quote"), header=None,
                         names=["ms", "bid", "bidsz", "ask", "asksz", "exch", "cond", "susp"])
    for df in (trades, quotes):
        df.query("ms >= @OPEN_MS and ms < @CLOSE_MS and susp == 0", inplace=True)
        df["sec"] = (df["ms"] - OPEN_MS) // 1000
    trades["price"] /= 1e4
    quotes["bid"] /= 1e4
    quotes["ask"] /= 1e4

    # NBBO at the end of each second from the latest quote of every exchange
    def best(side, agg):
        q = quotes.loc[quotes[side] > 0, ["sec", "exch", side]]
        last = q.groupby(["sec", "exch"])[side].last().unstack("exch")
        last = last.reindex(range(N_SEC)).ffill()
        return getattr(last, agg)(axis=1).to_numpy()

    bid, ask = best("bid", "max"), best("ask", "min")
    mid = 0.5 * (bid + ask)
    bad = ~np.isfinite(mid) | (bid >= ask)
    mid[bad] = np.nan
    mid = pd.Series(mid).ffill().bfill().to_numpy()
    spread_bp = 1e4 * (ask - bid) / mid
    spread_bp = pd.Series(np.where(bad, np.nan, spread_bp)).ffill().bfill().to_numpy()

    # sign trades against the midquote prevailing at the end of the previous second
    prev_mid = np.concatenate([[mid[0]], mid[:-1]])
    ref = prev_mid[trades["sec"].to_numpy()]
    price = trades["price"].to_numpy()
    sign = np.sign(price - ref)
    tick = np.sign(np.diff(price, prepend=price[0]))
    tick = pd.Series(tick).replace(0, np.nan).ffill().fillna(0).to_numpy()
    sign = np.where(sign == 0, tick, sign)
    sec = trades["sec"].to_numpy()
    size = trades["size"].to_numpy()
    out = {
        "mid": mid,
        "ret_bp": np.concatenate([[0.0], 1e4 * np.diff(np.log(mid))]),
        "spread_bp": spread_bp,
        "n_buy": np.bincount(sec, weights=(sign > 0), minlength=N_SEC),
        "n_sell": np.bincount(sec, weights=(sign < 0), minlength=N_SEC),
        "vol_buy": np.bincount(sec, weights=size * (sign > 0), minlength=N_SEC),
        "vol_sell": np.bincount(sec, weights=size * (sign < 0), minlength=N_SEC),
        "volume": np.bincount(sec, weights=size, minlength=N_SEC),
        "n_trades": np.bincount(sec, minlength=N_SEC).astype(float),
    }
    return out


def load_symbol(sym="spy", days=DAYS):
    """Per-second series for all days, concatenated, with a day index."""
    cache = os.path.join(DATA_DIR, f"micro_{sym}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    parts = [process_day(sym, d) for d in days]
    data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    data["day"] = np.repeat(np.arange(len(days)), N_SEC)
    data["sec"] = np.tile(np.arange(N_SEC), len(days))
    np.savez(cache, **data)
    return data


if __name__ == "__main__":
    import sys

    for sym in sys.argv[1:] or ["spy"]:
        d = load_symbol(sym)
        print(f"{sym}: {d['mid'].size} seconds over {len(DAYS)} days; mean trades/s {d['n_trades'].mean():.1f}, "
              f"mean volume/s {d['volume'].mean():.0f}, mean spread {d['spread_bp'].mean():.2f} bp, "
              f"return sd {d['ret_bp'].std():.2f} bp, imbalance sd {(d['n_buy'] - d['n_sell']).std():.2f}")
