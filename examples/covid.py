"""COVID-19 in the United States: daily deaths as a distributed lag of cases.

Deaths follow infections with a delay of a few weeks, so the lag curve relating
daily deaths to earlier case counts estimates both the size of the effect and
the shape of the delay distribution. The response is a count, so a negative
binomial model with a log link is used; weekday indicators absorb the strong
day-of-week pattern in death reporting. The exposure is the log of the
seven-day trailing mean of new cases.

Data: The New York Times, cumulative US cases and deaths; June 2020 to March
2021, a period with comparatively stable testing and no vaccination effect on
the case fatality ratio until its end.

Run: python examples/covid.py
"""
import datetime as dt
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import fetch, read_table  # noqa: E402

from badlm import BADLM  # noqa: E402
from badlm.plotting import plot_lag_curve  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
paperfig = os.path.join(here, "..", "paper", "figures")

# ---------------------------------------------------------------- data
raw = read_table(fetch("covid-us.csv"))
dates = [dt.date.fromisoformat(d) for d in raw["date"]]
cases = np.maximum(np.diff(np.array(raw["cases"])), 0.0)
deaths = np.maximum(np.diff(np.array(raw["deaths"])), 0.0)
dates = dates[1:]
start, end = dt.date(2020, 6, 1), dt.date(2021, 3, 31)
p = 35
i0 = dates.index(start) - p                      # keep p days of history for the lags
i1 = dates.index(end) + 1
dates, cases, deaths = dates[i0:i1], cases[i0:i1], deaths[i0:i1]
cases7 = np.convolve(cases, np.ones(7) / 7, mode="full")[: len(cases)]   # trailing mean
x = np.log1p(cases7)
y = deaths
weekday = np.array([d.weekday() for d in dates])
W = np.stack([(weekday == k).astype(float) for k in range(6)], axis=1)   # Mon..Sat; Sunday baseline
print(f"{len(y) - p} days of deaths from {dates[p]} to {dates[-1]}, mean {y[p:].mean():.0f} per day")

# ---------------------------------------------------------------- fit
t0 = time.time()
model = BADLM(nlag=p, k=24, family="negbin", end_constraint="zero")
model.fit(x, y, covariates=W, samples=1000, warmup=1000, chains=2, seed=3)
print(f"fitted in {time.time() - t0:.0f}s")
print(model.summary())
ce = model.cumulative_effect()
gamma = model.draws_["theta"][..., model.struct.K_r + 1:].reshape(-1, 6)
names = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]
print("weekday effects relative to Sunday (multiplicative):",
      ", ".join(f"{n} {np.exp(g):.2f}" for n, g in zip(names, gamma.mean(axis=0))))
print(f"dispersion (NB concentration) mean {model.draws_['nb_conc'].mean():.1f}")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(14, 4), gridspec_kw={"width_ratios": [1.1, 1.1, 1.3]})
t = np.array(dates[p:])
axes[0].plot(t, cases7[p:] / 1e3, color="C0", lw=1.5)
axes[0].set_title("(a) new cases, 7-day mean (thousands)")
axes[1].plot(t, y[p:], color="C0", lw=0.8)
axes[1].set_title("(b) daily deaths")
for ax in axes[:2]:
    ax.tick_params(axis="x", labelrotation=30)
plot_lag_curve(model.summary_, ax=axes[2], title="(c) elasticity of deaths to lagged cases")
axes[2].set_xlabel("lag (days)")
axes[2].set_ylabel("coefficient")
fig.tight_layout()
for d in (figdir, paperfig):
    os.makedirs(d, exist_ok=True)
fig.savefig(os.path.join(figdir, "covid_lag_curve.png"), dpi=130)
fig.savefig(os.path.join(paperfig, "covid_lag_curve.pdf"))
print("figures written")
