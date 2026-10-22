"""Daily market data for the trading experiments: a cross-asset panel of
adjusted closes (and volume where available) from public GitHub sources.

* US equities and ETFs 1998-2021: daily bars from the QuantConnect Lean sample
  data, adjusted for splits and dividends with Lean's factor files.
* 18 floating exchange rates against the US dollar 1990-2026 (Federal Reserve
  H.10 via the frictionlessdata ``exchange-rates`` package), expressed as the
  dollar value of one unit of foreign currency so that a positive return is a
  gain to a dollar investor holding the currency.
* Spot WTI and Brent (``oil-prices``), gold from the Lean XAUUSD CFD series, and the CBOE VIX (``finance-vix``) as a
  state variable.

Everything is cached in ``examples/data/markets.pkl``.
"""
from __future__ import annotations

import io
import os
import pickle
import urllib.request
import zipfile

import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
LEAN = "https://raw.githubusercontent.com/QuantConnect/Lean/master/Data/"
FD = "https://raw.githubusercontent.com/datasets/"
# Lean's QQQ file is empty between 2004 and March 2011, so it starts after the gap
EQUITIES = {"spy": "1998-01-01", "qqq": "2011-03-24", "iwm": "2000-05-26", "eem": "2003-04-14",
            "uso": "2006-04-10", "aapl": "1998-01-01", "ibm": "1998-01-01", "bac": "1998-01-01",
            "aig": "1998-01-01", "goog": "2014-04-04", "fb": "2012-05-18"}
SPLICES = {}
# the frictionlessdata package quotes every currency as foreign units per dollar
FX = {"Australia": False, "Canada": False, "Denmark": False, "Euro": False, "Japan": False,
      "New Zealand": False, "Norway": False, "Sweden": False, "Switzerland": False,
      "United Kingdom": False, "Brazil": False, "Mexico": False, "South Africa": False,
      "South Korea": False, "India": False, "Singapore": False, "Taiwan": False, "Thailand": False}
# managed or stale before 1995 (many days without a change)
EM_FROM_1995 = ("BRL", "MXN", "ZAR", "KRW", "INR", "TWD", "THB")
FX_CODE = {"Australia": "AUD", "Canada": "CAD", "Denmark": "DKK", "Euro": "EUR", "Japan": "JPY",
           "New Zealand": "NZD", "Norway": "NOK", "Sweden": "SEK", "Switzerland": "CHF",
           "United Kingdom": "GBP", "Brazil": "BRL", "Mexico": "MXN", "South Africa": "ZAR",
           "South Korea": "KRW", "India": "INR", "Singapore": "SGD", "Taiwan": "TWD", "Thailand": "THB"}


def _get(url, name):
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, name)
    if not os.path.exists(path):
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, path)
    return path


def _lean_zip_csv(path, names):
    with zipfile.ZipFile(path) as zf:
        raw = zf.read(zf.namelist()[0])
    df = pd.read_csv(io.BytesIO(raw), header=None, names=names)
    df["date"] = pd.to_datetime(df["date"].astype(str).str[:8], format="%Y%m%d")
    return df.set_index("date")


def lean_equity(sym):
    bars = _lean_zip_csv(_get(LEAN + f"equity/usa/daily/{sym}.zip", f"{sym}_daily.zip"),
                         ["date", "open", "high", "low", "close", "volume"])
    bars[["open", "high", "low", "close"]] /= 1e4
    fac = pd.read_csv(_get(LEAN + f"equity/usa/factor_files/{sym}.csv", f"{sym}_factors.csv"), header=None,
                      names=["date", "price_factor", "split_factor", "ref"], usecols=[0, 1, 2, 3])
    fac = fac[fac["date"].astype(str).str.len() == 8]
    fac["date"] = pd.to_datetime(fac["date"].astype(str), format="%Y%m%d")
    fac = fac.sort_values("date")
    # a bar on day d uses the factors of the first factor row dated on or after d
    idx = np.searchsorted(fac["date"].to_numpy(), bars.index.to_numpy(), side="left")
    idx = np.minimum(idx, len(fac) - 1)
    adj = fac["price_factor"].to_numpy()[idx] * fac["split_factor"].to_numpy()[idx]
    out = pd.DataFrame({"close": bars["close"].to_numpy() * adj,
                        "volume": bars["volume"].to_numpy() / fac["split_factor"].to_numpy()[idx]},
                       index=bars.index)
    out = out[~out.index.duplicated()].loc[EQUITIES[sym]:]
    if sym in SPLICES:
        d = pd.Timestamp(SPLICES[sym])
        before = out.index < d
        prev = out.loc[before, "close"].iloc[-1]
        out.loc[before, "close"] *= out.loc[d, "close"] / prev
    return out


def fred_fx():
    df = pd.read_csv(_get(FD + "exchange-rates/main/data/daily.csv", "fx_daily.csv"))
    df["Date"] = pd.to_datetime(df["Date"])
    out = {}
    for country, usd_per_unit in FX.items():
        s = df[df["Country"] == country].set_index("Date")["Exchange rate"]
        s = pd.to_numeric(s, errors="coerce").dropna()
        s = s[s > 0]
        first = "1995-01-01" if FX_CODE[country] in EM_FROM_1995 else "1990-01-01"
        out[FX_CODE[country]] = (s if usd_per_unit else 1.0 / s).loc[first:]
    return pd.DataFrame(out)


def commodities():
    out = {}
    # Henry Hub spot is stale on a third of days and jumps by over 100% on others: excluded
    for name, url in (("WTI", "oil-prices/main/data/wti-daily.csv"), ("BRENT", "oil-prices/main/data/brent-daily.csv")):
        d = pd.read_csv(_get(FD + url, os.path.basename(url).replace("daily", name.lower())))
        s = pd.Series(pd.to_numeric(d["Price"], errors="coerce").to_numpy(), index=pd.to_datetime(d["Date"]))
        out[name] = s[s > 0].dropna()
    gold = _lean_zip_csv(_get(LEAN + "cfd/oanda/daily/xauusd.zip", "xauusd_daily.zip"),
                         ["date", "bo", "bh", "bl", "bc", "bv", "ao", "ah", "al", "ac", "av"])
    out["GOLD"] = 0.5 * (gold["bc"] + gold["ac"])
    return pd.DataFrame(out)


def vix():
    d = pd.read_csv(_get(FD + "finance-vix/main/data/vix-daily.csv", "vix-daily.csv"))
    return pd.Series(d["CLOSE"].to_numpy(), index=pd.to_datetime(d["DATE"]), name="VIX")


def load_panel(refresh=False):
    """Dict with ``close`` (DataFrame, business days, forward-filled within each
    series' own history), ``volume`` (equities only), ``vix`` and ``kind``
    (series name -> 'equity' | 'fx' | 'commodity')."""
    os.makedirs(DATA_DIR, exist_ok=True)
    cache = os.path.join(DATA_DIR, "markets.pkl")
    if os.path.exists(cache) and not refresh:
        with open(cache, "rb") as f:
            return pickle.load(f)
    closes, volumes, kind = {}, {}, {}
    for sym in EQUITIES:
        eq = lean_equity(sym)
        closes[sym.upper()] = eq["close"]
        volumes[sym.upper()] = eq["volume"]
        kind[sym.upper()] = "equity"
    fx = fred_fx()
    for c in fx:
        closes[c] = fx[c]
        kind[c] = "fx"
    cm = commodities()
    for c in cm:
        closes[c] = cm[c]
        kind[c] = "commodity"
    close = pd.DataFrame(closes).sort_index()
    close = close[close.index.dayofweek < 5]
    # forward-fill holidays inside each series' own history only
    inside = close.notna().cummax() & close.notna()[::-1].cummax()[::-1]
    close = close.ffill().where(inside)
    volume = pd.DataFrame(volumes).reindex(close.index)
    panel = {"close": close, "volume": volume, "vix": vix().reindex(close.index).ffill(), "kind": kind}
    with open(cache, "wb") as f:
        pickle.dump(panel, f)
    return panel


def log_returns(close):
    return np.log(close).diff()


def ewma_vol(ret, halflife=30, min_periods=60, floor_frac=0.25):
    """Ex-ante daily volatility: exponentially weighted standard deviation of
    past returns, lagged one day so that it is known when the position is set,
    floored at ``floor_frac`` of the series' expanding median volatility so
    that a quiet spell cannot produce an enormous position."""
    v = ret.ewm(halflife=halflife, min_periods=min_periods).std().shift(1)
    floor = floor_frac * v.expanding().median()
    return v.clip(lower=floor)


if __name__ == "__main__":
    p = load_panel(refresh=True)
    c = p["close"]
    r = log_returns(c)
    print(f"{c.shape[1]} series, {c.index[0].date()} to {c.index[-1].date()}")
    print("| Series | Kind | First | Last | Days | Ann. vol % | Ann. mean % |")
    print("|---|---|---|---|---|---|---|")
    for s in c:
        rs = r[s].dropna()
        print(f"| {s} | {p['kind'][s]} | {rs.index[0].date()} | {rs.index[-1].date()} | {rs.size} | "
              f"{100*rs.std()*np.sqrt(252):.1f} | {100*rs.mean()*252:.1f} |")
