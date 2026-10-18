"""Lag curves for daily stock market data.

Prices are transformed to percentage log returns, ``r_t = 100 log(P_t / P_{t-1})``,
and four questions are asked of the S&P 500, Apple and the VIX, 2007-2016:

(a) does the next day's S&P 500 return depend on the last 20 days of returns?
(b) the same for Apple;
(c) how does implied volatility (the daily change in log VIX) respond to S&P
    500 returns on the same day and the 60 days before;
(d) how does today's log absolute return (a measure of realised volatility)
    depend on the last 100 days of log absolute returns (volatility clustering).

Run: python examples/stocks.py
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import fetch, read_table  # noqa: E402

from badlm import BADLM, LagCurveSummary  # noqa: E402
from badlm.plotting import plot_lag_curve  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
paperfig = os.path.join(here, "..", "paper", "figures")

# ---------------------------------------------------------------- data
stocks = read_table(fetch("stockdata.csv"))
vix = read_table(fetch("vix-daily.csv"))
vix_close = dict(zip(vix["DATE"], vix["CLOSE"]))
keep = [i for i, d in enumerate(stocks["Date"]) if d in vix_close]
dates = [stocks["Date"][i] for i in keep]
sp = np.array([stocks["GSPC"][i] for i in keep])
aapl = np.array([stocks["AAPL"][i] for i in keep])
vx = np.array([vix_close[d] for d in dates])

r_sp = 100.0 * np.diff(np.log(sp))        # percentage log returns, days 1..T
r_aapl = 100.0 * np.diff(np.log(aapl))
d_lvix = 100.0 * np.diff(np.log(vx))      # percentage change in VIX, aligned with r_sp
print(f"{len(dates)} trading days from {dates[0]} to {dates[-1]}; "
      f"S&P 500 daily return sd {r_sp.std():.2f}%, Apple {r_aapl.std():.2f}%, "
      f"VIX change sd {d_lvix.std():.1f}%")

# ---------------------------------------------------------------- fits
mcmc = dict(samples=1000, warmup=1000, chains=2, seed=1)
t0 = time.time()
fits = {}
# (a), (b): y_t = r_{t+1} on x_t = r_t, so lag j of x is lag j + 1 of the return
fits["a"] = BADLM(nlag=20, k=14).fit(r_sp[:-1], r_sp[1:], **mcmc)
fits["b"] = BADLM(nlag=20, k=14).fit(r_aapl[:-1], r_aapl[1:], **mcmc)
# (c): contemporaneous and lagged response of implied volatility to returns
fits["c"] = BADLM(nlag=60, k=40).fit(r_sp, d_lvix, **mcmc)
# (d): volatility clustering, log|r_{t+1}| on log|r_t|, ..., log|r_{t-100}|; the log
# transform tames the heavy tails of absolute returns (0.1% added to avoid log 0)
la_sp = np.log(np.abs(r_sp) + 0.1)
fits["d"] = BADLM(nlag=100, k=68).fit(la_sp[:-1], la_sp[1:], **mcmc)
print(f"four fits in {time.time() - t0:.0f}s")

titles = {
    "a": "(a) S&P 500: next-day return on past returns",
    "b": "(b) Apple: next-day return on past returns",
    "c": "(c) Change in log VIX on S&P 500 returns",
    "d": "(d) S&P 500: log absolute return on its own past",
}
xlab = {"a": "lag (days before the predicted day)", "b": "lag (days before the predicted day)",
        "c": "lag (days)", "d": "lag (days before the predicted day)"}

for key, lag in (("a", 0), ("c", 0), ("c", 1), ("d", 0)):
    s = fits[key].summary_
    print(f"panel {key}, x lag {lag}: median {s.median[lag]:.3f} ({s.lower[lag]:.3f}, {s.upper[lag]:.3f})")
print("\n| Panel | Cumulative effect (95% interval) | ED | Peak lag | Effective lag length | max R-hat | Divergences |")
print("|---|---|---|---|---|---|---|")
for key, m in fits.items():
    ce = m.cumulative_effect()
    d = m.diagnostics_
    rh = max(np.nanmax(v["rhat"]) for v in d.values() if isinstance(v, dict))
    print(f"| {key} | {ce['mean']:.3f} ({ce['lower']:.3f}, {ce['upper']:.3f}) | {m.effective_df()['mean']:.1f} | "
          f"{m.peak_lag()['mode']} | {m.effective_lag_length()['length']} | {rh:.2f} | {d.get('divergences', 0)} |")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(2, 2, figsize=(11, 7))
for ax, key in zip(axes.ravel(), fits):
    m = fits[key]
    s = m.summary_
    # in (a), (b), (d) the exposure is the previous day's return, so lag j of x
    # is lag j + 1 of the return relative to the predicted day
    shift = 1 if key in ("a", "b", "d") else 0
    shifted = LagCurveSummary(lag=s.lag + shift, median=s.median, lower=s.lower, upper=s.upper,
                              mean=s.mean, probs=s.probs)
    plot_lag_curve(shifted, ax=ax, title=titles[key])
    ax.set_xticks(np.linspace(shifted.lag[0], shifted.lag[-1], 6).round())
    ax.set_xlabel(xlab[key])
    ax.set_ylabel("elasticity" if key == "d" else "coefficient (% per %)")
    ax.get_legend().remove()
axes[0, 0].legend(frameon=False, loc="upper right")
fig.tight_layout()
for d in (figdir, paperfig):
    os.makedirs(d, exist_ok=True)
fig.savefig(os.path.join(figdir, "stocks_lag_curves.png"), dpi=130)
fig.savefig(os.path.join(paperfig, "stocks_lag_curves.pdf"))
print("figures written")
