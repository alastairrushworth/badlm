"""Plotting helpers."""
from __future__ import annotations

import numpy as np


def plot_lag_curve(summary, ax=None, truth=None, color="C0", label="posterior median",
                   title="Posterior median and 95% credible band for the lag curve"):
    """Ribbon plot of the lag-curve summary returned by ``BADLM.lag_curve_summary``."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    ax.fill_between(summary.lag, summary.lower, summary.upper, color=color, alpha=0.25, lw=0,
                    label="95% credible band")
    ax.plot(summary.lag, summary.median, color=color, lw=2, label=label)
    if truth is not None:
        ax.plot(summary.lag, np.asarray(truth), color="k", lw=1, ls="--", label="truth")
    ax.axhline(0.0, color="grey", lw=0.5)
    ax.set_xlabel("lag")
    ax.set_ylabel("lag coefficient")
    ax.set_title(title)
    ax.legend(frameon=False)
    return ax
