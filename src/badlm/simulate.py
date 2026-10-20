"""Simulated distributed-lag data and reference lag curves."""
from __future__ import annotations

import numpy as np

from .design import lag_matrix


def ar1(n, phi, sd, rng, burn=100):
    e = rng.normal(scale=sd, size=n + burn)
    x = np.empty(n + burn)
    x[0] = e[0] / np.sqrt(1.0 - phi**2)
    for t in range(1, n + burn):
        x[t] = phi * x[t - 1] + e[t]
    return x[burn:]


def simulate_dlm(lag_curve, n=500, ar_x=0.5, sd_x=0.1, ar_eps=0.2, sd_eps=0.1,
                 standardize_x=False, seed=None):
    """Simulate an exposure series and a Gaussian response under a lag curve.

    Follows the design of the paper's simulation study: ``x_t ~ N(ar_x x_{t-1},
    sd_x^2)`` and ``y_t = sum_j beta_j x_{t-j} + eps_t`` with AR(1) errors.
    Returns ``(x, y)`` where ``y`` has length ``n - p``.
    """
    rng = np.random.default_rng(seed)
    beta = np.asarray(lag_curve, float)
    p = beta.size - 1
    x = ar1(n, ar_x, sd_x, rng)
    if standardize_x:
        x = (x - x.mean()) / x.std(ddof=1)
    eps = ar1(n - p, ar_eps, sd_eps, rng)
    y = lag_matrix(x, p) @ beta + eps
    return x, y


def readme_lag_curve(p=50):
    """The 'hump' lag function from the original package README."""
    lag = np.arange(p + 1, dtype=float)
    pdf = np.exp(-0.5 * ((lag - 10.0) / 4.0) ** 2) / (4.0 * np.sqrt(2 * np.pi))
    return -0.1 + (0.01 * np.exp(-0.2 * lag) + np.exp(pdf)) / 10.0


def paper_lag_curves(p=50):
    """Five lag-curve shapes in the spirit of Figure 1 of the paper.

    The published curves were not released in closed form; these are
    approximations with the same qualitative features and scale.
    """
    lag = np.arange(p + 1, dtype=float)
    bump = lambda c, s: np.exp(-0.5 * ((lag - c) / s) ** 2)  # noqa: E731
    return {
        "delayed_peak": 0.022 * bump(11.0, 2.2) + 0.001 * np.exp(-lag / 4.0),
        "decay_curve": 0.02 * np.exp(-lag / 8.0),
        "flat_response": np.zeros(p + 1),
        "displacement": 0.018 * bump(9.0, 2.5) - 0.006 * bump(16.0, 4.0) * (lag > 10),
        "sharp_peak": 0.04 * np.exp(-lag / 1.2),
    }
