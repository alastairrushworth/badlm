"""The wide panel for the multivariate experiments, and the frozen held-out
period.

Sources (all free, fetched on first use and cached in ``examples/data/``):

* ``skfolio`` (PyPI wheel, bundled data): adjusted daily closes of 20 S&P 500
  constituents and the index, 1990-2022; five factor ETFs 2014-2022.
* ``skfolio/skfolio-datasets`` (GitHub): 64 FTSE 100 constituents 2000-2023
  in pounds, 1,455 Nasdaq names 2018-2023, and 3-month implied volatilities
  of the 20 S&P 500 names 2010-2022.
* ``robertmartin8/PyPortfolioOpt`` test resources: 20 US stocks 1990-2018.
* The currencies, oil, gold and Lean equities of ``_markets.py``.

``HOLDOUT_START`` is the first day of the held-out period. Development runs
stop before it; it is used once, at the end, by passing ``--final``.
"""
from __future__ import annotations

import glob
import io
import os
import subprocess
import sys
import urllib.request
import zipfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _markets import DATA_DIR, load_panel  # noqa: E402

HOLDOUT_START = pd.Timestamp("2019-01-01")
SECTOR = {  # GICS-style sector of the US names, for industry baskets
    "AAPL": "tech", "AMD": "tech", "MSFT": "tech", "GOOG": "tech", "FB": "tech", "IBM": "tech", "AMZN": "cons",
    "BAC": "fin", "JPM": "fin", "MA": "fin", "AIG": "fin", "BBY": "cons", "HD": "cons", "WMT": "cons",
    "SBUX": "cons", "GM": "cons", "UAA": "cons", "SHLD": "cons", "CVX": "energy", "XOM": "energy",
    "RRC": "energy", "GE": "ind", "LLY": "health", "MRK": "health", "PFE": "health", "JNJ": "health",
    "KO": "cons", "PEP": "cons", "PG": "cons", "T": "tel", "UNH": "health", "WMT": "cons", "XOM": "energy",
    "SPY": "index", "QQQ": "index", "IWM": "index", "EEM": "index", "USO": "energy",
}


def _skfolio_bundled():
    """The three datasets inside the skfolio wheel."""
    want = ("sp500_dataset", "sp500_index", "factors_dataset")
    if all(os.path.exists(os.path.join(DATA_DIR, f"{w}.csv.gz")) for w in want):
        return
    tmp = os.path.join(DATA_DIR, "_skfolio_wheel")
    os.makedirs(tmp, exist_ok=True)
    subprocess.run([sys.executable, "-m", "pip", "download", "-q", "--no-deps", "skfolio", "-d", tmp], check=True)
    whl = glob.glob(os.path.join(tmp, "*.whl"))[0]
    with zipfile.ZipFile(whl) as zf:
        for name in zf.namelist():
            base = os.path.basename(name)
            if base.replace(".csv.gz", "") in want:
                with open(os.path.join(DATA_DIR, base), "wb") as f:
                    f.write(zf.read(name))


def _skfolio_remote():
    want = ("ftse100_dataset", "nasdaq_dataset", "sp500_implied_vol_dataset")
    if all(os.path.exists(os.path.join(DATA_DIR, f"{w}.csv.gz")) for w in want):
        return
    tmp = os.path.join(DATA_DIR, "_skfolio_datasets")
    subprocess.run(["git", "clone", "-q", "--depth", "1", "https://github.com/skfolio/skfolio-datasets", tmp], check=True)
    for w in want:
        os.replace(os.path.join(tmp, "datasets", f"{w}.csv.gz"), os.path.join(DATA_DIR, f"{w}.csv.gz"))


def _read(name):
    return pd.read_csv(os.path.join(DATA_DIR, name), index_col=0, parse_dates=True)


def load_wide_panel():
    """Dict with ``close`` (business days), ``kind`` (equity_us, equity_uk,
    index, factor_etf, fx, commodity) and ``currency`` per series."""
    os.makedirs(DATA_DIR, exist_ok=True)
    _skfolio_bundled()
    _skfolio_remote()
    ppo = os.path.join(DATA_DIR, "ppo_stock_prices.csv")
    if not os.path.exists(ppo):
        urllib.request.urlretrieve(
            "https://raw.githubusercontent.com/robertmartin8/PyPortfolioOpt/master/tests/resources/stock_prices.csv", ppo)
    closes, kind, cur = {}, {}, {}
    sp = _read("sp500_dataset.csv.gz")
    for c in sp:
        closes[c] = sp[c]; kind[c] = "equity_us"; cur[c] = "USD"
    pp = pd.read_csv(ppo, index_col=0, parse_dates=True)
    for c in pp:
        if c not in closes and c not in ("SHLD",):
            closes[c] = pp[c]; kind[c] = "equity_us"; cur[c] = "USD"
    base = load_panel()
    for c in base["close"]:
        if c in closes:
            continue
        k = base["kind"][c]
        closes[c] = base["close"][c]
        kind[c] = {"equity": "index" if c in ("SPY", "QQQ", "IWM", "EEM") else "equity_us", "fx": "fx", "commodity": "commodity"}[k]
        cur[c] = "USD"
    idx = _read("sp500_index.csv.gz")
    closes["SP500"] = idx["SP500"]; kind["SP500"] = "index"; cur["SP500"] = "USD"
    fe = _read("factors_dataset.csv.gz")
    for c in fe:
        closes[c] = fe[c]; kind[c] = "factor_etf"; cur[c] = "USD"
    ft = _read("ftse100_dataset.csv.gz")
    for c in ft:
        closes[c] = ft[c]; kind[c] = "equity_uk"; cur[c] = "GBP"
    close = pd.DataFrame(closes).sort_index()
    close = close[close.index.dayofweek < 5]
    inside = close.notna().cummax() & close.notna()[::-1].cummax()[::-1]
    close = close.ffill().where(inside)
    return {"close": close, "kind": kind, "currency": cur}


def load_nasdaq():
    _skfolio_remote()
    return _read("nasdaq_dataset.csv.gz")


def load_implied_vol():
    _skfolio_remote()
    return _read("sp500_implied_vol_dataset.csv.gz")


if __name__ == "__main__":
    p = load_wide_panel()
    c = p["close"]
    r = np.log(c).diff()
    print(f"{c.shape[1]} series, {c.index[0].date()} to {c.index[-1].date()}; held out from {HOLDOUT_START.date()}")
    by_kind = pd.Series(p["kind"]).value_counts()
    print(by_kind.to_string())
    stale = (r == 0).mean()
    print("series with over 10% zero-return days:", list(stale[stale > 0.10].index))
    big = r.abs().max()
    print("series with a daily move over 50%:", list(big[big > 0.5].index))
    n = load_nasdaq(); iv = load_implied_vol()
    print("nasdaq", n.shape, "implied vol", iv.shape, iv.iloc[-1, :3].round(3).to_dict())
