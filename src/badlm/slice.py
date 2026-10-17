"""Univariate slice sampler with stepping out and shrinkage (Neal, 2003)."""
from __future__ import annotations

import numpy as np


def slice_sample(logf, x0, rng, width=1.0, max_steps=50, logf0=None):
    """One update of a scalar variable with log density ``logf``.

    Returns ``(x_new, logf(x_new))``. ``logf0`` may be passed to avoid
    recomputing the density at ``x0``.
    """
    if logf0 is None:
        logf0 = logf(x0)
    logy = logf0 + np.log(rng.uniform())
    u = rng.uniform()
    left = x0 - width * u
    right = left + width
    j = int(np.floor(max_steps * rng.uniform()))
    k = max_steps - 1 - j
    while j > 0 and logy < logf(left):
        left -= width
        j -= 1
    while k > 0 and logy < logf(right):
        right += width
        k -= 1
    while True:
        x1 = rng.uniform(left, right)
        lf1 = logf(x1)
        if logy < lf1:
            return x1, lf1
        if x1 < x0:
            left = x1
        else:
            right = x1
