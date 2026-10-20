"""Download and parse the public datasets used by the examples.

Files are cached in ``examples/data/`` (not committed). Sources:

* ``stockdata.csv``: daily closes of MSFT, IBM, SBUX, AAPL and the S&P 500
  (GSPC), 2007-2016, from the plotly/datasets repository.
* ``vix-daily.csv``: daily CBOE VIX from the datasets/finance-vix repository
  (Open Data Commons PDDL).
* ``covid-us.csv``: cumulative US COVID-19 cases and deaths compiled by The New
  York Times (nytimes/covid-19-data), used here for non-commercial illustration.
"""
from __future__ import annotations

import csv
import os
import urllib.request

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
SOURCES = {
    "stockdata.csv": "https://raw.githubusercontent.com/plotly/datasets/master/stockdata.csv",
    "vix-daily.csv": "https://raw.githubusercontent.com/datasets/finance-vix/main/data/vix-daily.csv",
    "covid-us.csv": "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us.csv",
}


def fetch(name):
    """Return the local path of a dataset, downloading it on first use."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, name)
    if not os.path.exists(path):
        print(f"downloading {SOURCES[name]}")
        urllib.request.urlretrieve(SOURCES[name], path)
    return path


def read_table(path):
    """Read a CSV into a dict of columns; numeric columns become floats."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]
    cols = {h: [r[i] for r in body] for i, h in enumerate(header)}
    for h, vals in cols.items():
        try:
            cols[h] = [float(v) for v in vals]
        except ValueError:
            pass
    return cols
