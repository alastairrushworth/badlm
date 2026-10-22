"""badlm: Bayesian adaptive distributed lag models.

The lag curve relating an exposure series to a response is represented by a
B-spline basis whose coefficients receive a first-difference penalty with
lag-varying precision; the log precisions are themselves smoothed, so the amount
of smoothing adapts to the shape of the curve (Rushworth, 2018,
arXiv:1801.06670).

The primary sampler integrates the spline coefficients out analytically and runs
the No-U-Turn Sampler on the remaining low-dimensional posterior; a Gibbs sampler
with slice updates and a Laplace approximation are also provided.
"""
from __future__ import annotations

import jax

# Collapsed likelihoods involve Cholesky factors of matrices whose condition
# number is governed by exp(tau); float64 is required for them to be reliable.
jax.config.update("jax_enable_x64", True)

from .design import bspline_basis, lag_matrix  # noqa: E402
from .model import BADLM, LagCurveSummary  # noqa: E402
from .simulate import paper_lag_curves, readme_lag_curve, simulate_dlm  # noqa: E402

__all__ = [
    "BADLM",
    "LagCurveSummary",
    "bspline_basis",
    "lag_matrix",
    "simulate_dlm",
    "paper_lag_curves",
    "readme_lag_curve",
]
__version__ = "0.3.0"
