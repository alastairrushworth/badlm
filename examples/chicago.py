"""Air pollution and mortality in Chicago, 1987 to 2000.

The classic distributed-lag setting of environmental epidemiology: daily deaths
regressed on particulate pollution (PM10) over the preceding weeks, adjusting
for long-term trend and season, day of week, temperature and dew point. The lag
curve shows whether the effect of a polluted day is immediate, delayed, and
whether it is followed by negative coefficients, the "mortality displacement"
or harvesting pattern in which deaths brought forward by pollution are missing
from later days.

Data: ``chicagoNMMAPS`` from the R package dlnm (Gasparrini), itself drawn from
the National Morbidity, Mortality and Air Pollution Study; 5,114 days of
all-cause, cardiovascular and respiratory deaths, PM10 (micrograms per cubic
metre), ozone, temperature and dew point. Missing PM10 values (5%) are
interpolated linearly.

Model: Gaussian on log deaths with AR(1) errors, 40 lags of PM10, a B-spline of
time with 7 degrees of freedom per year, weekday indicators, and B-splines of
temperature and dew point, all unpenalised. Coefficients are reported as the
percentage change in deaths per 10 micrograms per cubic metre.

Run: python examples/chicago.py
"""
import csv
import gzip
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from badlm import BADLM, bspline_basis  # noqa: E402
from badlm.plotting import plot_lag_curve  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
paperfig = os.path.join(here, "..", "paper", "figures")

# ---------------------------------------------------------------- data
with gzip.open(os.path.join(here, "data", "chicagoNMMAPS.csv.gz"), "rt") as f:
    rows = list(csv.DictReader(f))
col = lambda k, cast=float: np.array([cast(r[k]) if r[k] != "NA" else np.nan for r in rows])  # noqa: E731
death, cvd, resp = col("death"), col("cvd"), col("resp")
temp, dptp, pm10 = col("temp"), col("dptp"), col("pm10")
dow = np.array([r["dow"] for r in rows])
n = death.size
idx = np.arange(n)
pm10 = np.interp(idx, idx[~np.isnan(pm10)], pm10[~np.isnan(pm10)])
print(f"{n} days, mean deaths {death.mean():.1f} (cardiovascular {cvd.mean():.1f}, respiratory {resp.mean():.1f}), "
      f"mean PM10 {pm10.mean():.1f}")

# confounders, unpenalised: trend and season, weekday, temperature, dew point
years = n / 365.25
W = np.hstack([
    bspline_basis(idx, nseg=int(round(7 * years)) - 3, deg=3, domain=(0, n - 1))[:, 1:],
    np.stack([(dow == d).astype(float) for d in ("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")], axis=1),
    bspline_basis(temp, nseg=3, deg=3)[:, 1:],
    bspline_basis(dptp, nseg=2, deg=3)[:, 1:],
])
print(f"{W.shape[1]} unpenalised covariate columns")

# ---------------------------------------------------------------- fits
p, k = 40, 28
scale = 1000.0          # log scale -> % change per 10 units
fits = {}
t0 = time.time()
for name, y in (("all-cause", death), ("cardiovascular", cvd), ("respiratory", resp)):
    m = BADLM(nlag=p, k=k, ar1=True).fit(pm10, np.log(y + 0.5), covariates=W,
                                          samples=1000, warmup=1000, chains=2, seed=1)
    fits[name] = m
    ce = m.cumulative_effect()
    d = m.diagnostics_
    rh = max(np.nanmax(v["rhat"]) for v in d.values() if isinstance(v, dict))
    s = m.summary_
    print(f"{name:15s} lag0 {scale*s.median[0]:+.2f}% ({scale*s.lower[0]:+.2f}, {scale*s.upper[0]:+.2f}); "
          f"lag1 {scale*s.median[1]:+.2f}%; cumulative {scale*ce['mean']:+.2f}% ({scale*ce['lower']:+.2f}, {scale*ce['upper']:+.2f}) "
          f"per 10 ug/m3; ED {m.effective_df()['mean']:.1f}; phi {m.draws_['phi'].mean():.2f}; max R-hat {rh:.2f}; "
          f"divergences {d.get('divergences', 0)}", flush=True)
print(f"three fits in {time.time() - t0:.0f}s")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(2, 2, figsize=(11, 7))
for ax, (name, m) in zip(axes.ravel()[:3], fits.items()):
    s = m.summary_
    from badlm import LagCurveSummary
    pct = LagCurveSummary(lag=s.lag, median=scale * s.median, lower=scale * s.lower, upper=scale * s.upper,
                          mean=scale * s.mean, probs=s.probs)
    plot_lag_curve(pct, ax=ax, title=f"{name} deaths")
    ax.set_xlabel("lag (days)")
    ax.set_ylabel("% change per 10 µg/m³ PM10")
    if name != "all-cause":
        ax.get_legend().remove()
ax = axes[1, 1]
for (name, m), color in zip(fits.items(), ("C0", "C1", "C2")):
    cum = scale * np.cumsum(m.lag_curve_.reshape(-1, p + 1), axis=1)
    q = np.quantile(cum, [0.025, 0.5, 0.975], axis=0)
    ax.plot(m.lags, q[1], color=color, lw=2, label=name)
    ax.fill_between(m.lags, q[0], q[2], color=color, alpha=0.15, lw=0)
ax.axhline(0, color="grey", lw=0.5)
ax.set_title("cumulative effect up to each lag")
ax.set_xlabel("lag (days)")
ax.set_ylabel("% change per 10 µg/m³ PM10")
ax.legend(frameon=False)
for a, lab in zip(axes.ravel(), "abcd"):
    a.set_title(f"({lab}) " + a.get_title())
fig.tight_layout()
for d in (figdir, paperfig):
    os.makedirs(d, exist_ok=True)
fig.savefig(os.path.join(figdir, "chicago_lag_curves.png"), dpi=130)
fig.savefig(os.path.join(paperfig, "chicago_lag_curves.pdf"))
print("figures written")
