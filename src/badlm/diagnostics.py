"""Convergence diagnostics for arrays of draws shaped ``(chains, draws, ...)``."""
from __future__ import annotations

import numpy as np
from numpyro.diagnostics import effective_sample_size, split_gelman_rubin


def rhat(draws):
    draws = np.asarray(draws)
    if draws.shape[0] < 2:
        return np.full(draws.shape[2:], np.nan)
    return np.asarray(split_gelman_rubin(draws))


def ess(draws):
    return np.asarray(effective_sample_size(np.asarray(draws)))


def summarise_draws(draws_dict, names=None):
    """R-hat and bulk effective sample size for every site in ``draws_dict``."""
    out = {}
    for name, arr in draws_dict.items():
        if names is not None and name not in names:
            continue
        arr = np.asarray(arr)
        if arr.ndim < 2 or arr.shape[1] < 4:
            continue
        out[name] = {"rhat": rhat(arr), "ess": ess(arr)}
    return out
