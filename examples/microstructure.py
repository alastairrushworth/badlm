"""Intraday microstructure of SPY at one-second resolution.

Five days of every trade and quote (October 2013) are reduced to per-second
series by ``_microstructure.py``: midquote returns in basis points, order-flow
imbalance (buyer-initiated minus seller-initiated trades), volume and the
quoted spread. Lags are formed within each trading day, and B-splines of time
of day absorb the intraday pattern in volume and volatility. Five lag curves
are fitted:

(a) price impact: midquote return on order-flow imbalance at lags 0 to 120 s;
    the cumulative sum of these coefficients, panel (b), is the response of the
    price to one unit of order flow as a function of time since the trade;
(c) memory of order flow: imbalance on its own past over 300 s;
(d) volatility on past volume: absolute return on volume over 120 s;
(e) volume on past volatility: volume on absolute returns over 120 s;
(f) memory of order flow for BAC, a large-tick stock, for comparison.

Run: python examples/microstructure.py   (requires pandas for the first run)
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _microstructure import N_SEC, load_symbol  # noqa: E402

from badlm import BADLM, LagCurveSummary, bspline_basis  # noqa: E402
from badlm.plotting import plot_lag_curve  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
figdir = os.path.join(here, "..", "figures")
paperfig = os.path.join(here, "..", "paper", "figures")


def shifted(x, y, day):
    """Within each day, pair x_t with y_{t+1}: the lag-j coefficient of x then
    refers to lag j + 1 of the predicted series."""
    keep_x = np.ones(x.size, bool)
    keep_x[np.where(np.diff(day) != 0)[0]] = False     # last second of each day
    keep_x[-1] = False
    keep_y = np.roll(keep_x, 1)
    return x[keep_x], y[keep_y], day[keep_x]


def tod_basis(sec, ncol=10):
    return bspline_basis(sec, nseg=ncol - 3, deg=3, domain=(0, N_SEC - 1))[:, 1:]


mcmc = dict(samples=1000, warmup=1000, chains=2, seed=1)
fits, info = {}, {}
t0 = time.time()
for sym in ("spy", "bac"):
    d = load_symbol(sym)
    day, sec = d["day"], d["sec"]
    imb = d["n_buy"] - d["n_sell"]
    ret = d["ret_bp"]
    vol = d["volume"] / 1e3
    tod = tod_basis(sec)
    print(f"{sym}: {ret.size} seconds; return sd {ret.std():.2f} bp; imbalance sd {imb.std():.2f} trades; "
          f"volume mean {vol.mean():.1f}k shares/s; spread {d['spread_bp'].mean():.2f} bp", flush=True)
    if sym == "spy":
        fits["a"] = BADLM(nlag=120, k=60).fit(imb, ret, covariates=tod, groups=day, **mcmc)
        xs, ys, gs = shifted(imb, imb, day)
        fits["c"] = BADLM(nlag=300, k=100).fit(xs, ys, covariates=tod[:xs.size], groups=gs, **mcmc)
        xs, ys, gs = shifted(vol, np.abs(ret), day)
        fits["d"] = BADLM(nlag=120, k=60).fit(xs, ys, covariates=tod[:xs.size], groups=gs, **mcmc)
        xs, ys, gs = shifted(np.abs(ret), vol, day)
        fits["e"] = BADLM(nlag=120, k=60).fit(xs, ys, covariates=tod[:xs.size], groups=gs, **mcmc)
    else:
        xs, ys, gs = shifted(imb, imb, day)
        fits["f"] = BADLM(nlag=300, k=100).fit(xs, ys, covariates=tod[:xs.size], groups=gs, **mcmc)
print(f"five fits in {time.time() - t0:.0f}s", flush=True)

titles = {"a": "(a) SPY: midquote return on order-flow imbalance",
          "b": "(b) SPY: cumulative price response to one unit of order flow",
          "c": "(c) SPY: order-flow imbalance on its own past",
          "d": "(d) SPY: absolute return on past volume",
          "e": "(e) SPY: volume on past absolute returns",
          "f": "(f) BAC: order-flow imbalance on its own past"}
units = {"a": "bp per trade of imbalance", "b": "bp per trade of imbalance", "c": "trades per trade",
         "d": "bp per 1000 shares", "e": "1000 shares per bp", "f": "trades per trade"}

print("\n| Panel | First-lag coefficient | Cumulative effect | ED | Effective lag length | max R-hat | Divergences |")
print("|---|---|---|---|---|---|---|")
for key in ("a", "c", "d", "e", "f"):
    m = fits[key]
    s = m.summary_
    ce = m.cumulative_effect()
    dg = m.diagnostics_
    rh = max(np.nanmax(v["rhat"]) for v in dg.values() if isinstance(v, dict))
    print(f"| {key} | {s.median[0]:.4f} ({s.lower[0]:.4f}, {s.upper[0]:.4f}) | {ce['mean']:.4f} ({ce['lower']:.4f}, {ce['upper']:.4f}) | "
          f"{m.effective_df()['mean']:.1f} | {m.effective_lag_length()['length']} | {rh:.2f} | {dg.get('divergences', 0)} |")
cum = np.cumsum(fits["a"].lag_curve_.reshape(-1, 121), axis=1)
q = np.quantile(cum, [0.025, 0.5, 0.975], axis=0)
for ell in (0, 1, 5, 30, 120):
    print(f"cumulative price response at {ell:3d} s: {q[1][ell]:.4f} ({q[0][ell]:.4f}, {q[2][ell]:.4f}) bp per trade")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
for ax, key in zip(axes.ravel(), "abcdef"):
    if key == "b":
        lags = np.arange(121)
        ax.fill_between(lags, q[0], q[2], color="C0", alpha=0.25, lw=0, label="95% credible band")
        ax.plot(lags, q[1], color="C0", lw=2, label="posterior median")
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_title(titles[key])
        ax.set_xlabel("seconds since the trade")
    else:
        m = fits[key]
        s = m.summary_
        shift = 0 if key == "a" else 1
        ss = LagCurveSummary(lag=s.lag + shift, median=s.median, lower=s.lower, upper=s.upper, mean=s.mean, probs=s.probs)
        plot_lag_curve(ss, ax=ax, title=titles[key])
        ax.set_xlabel("lag (seconds)" if key == "a" else "lag (seconds before the predicted second)")
        ax.get_legend().remove()
    ax.set_ylabel(units[key])
axes[0, 1].legend(frameon=False)
fig.tight_layout()
for dd in (figdir, paperfig):
    os.makedirs(dd, exist_ok=True)
fig.savefig(os.path.join(figdir, "microstructure_lag_curves.png"), dpi=130)
fig.savefig(os.path.join(paperfig, "microstructure_lag_curves.pdf"))
print("figures written")
