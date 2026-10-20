"""Prior densities shared by the NumPy Gibbs sampler and the tests."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Priors:
    """Hyperparameters of the model.

    * ``sigma2 ~ InverseGamma(sigma2_a, sigma2_b)``
    * log precisions: the level ``mean(tau)`` has a penalised-complexity prior,
      ``exp(-mean(tau) / 2) ~ Exponential(tau_level_rate)`` with
      ``tau_level_rate = -log(tau_level_tail_prob) / tau_level_scale``, whose
      density decays like ``exp(-tau / 2)`` for large precisions and
      double-exponentially for small ones; first differences
      ``tau_{k+1} - tau_k ~ N(0, zeta^2)``
    * ``zeta ~ Exponential(zeta_rate)``, a penalised-complexity prior that shrinks
      towards the non-adaptive model; ``zeta_rate = -log(zeta_tail_prob) / zeta_scale``
      gives ``P(zeta > zeta_scale) = zeta_tail_prob``.
    * ``tau_basis`` number of B-spline functions for ``tau_model='basis'``
    * ``hs_scale`` global scale of the horseshoe alternative
    * ``coef_sd`` prior sd of unpenalised coefficients in non-Gaussian models
    """

    sigma2_a: float = 0.001
    sigma2_b: float = 0.001
    tau_level_scale: float = 1.0
    tau_level_tail_prob: float = 0.05
    zeta_scale: float = 2.0
    zeta_tail_prob: float = 0.1
    tau_basis: int = 6
    hs_scale: float = 1.0
    coef_sd: float = 10.0

    @property
    def zeta_rate(self):
        return -np.log(self.zeta_tail_prob) / self.zeta_scale

    @property
    def tau_level_rate(self):
        return -np.log(self.tau_level_tail_prob) / self.tau_level_scale

    @property
    def tau_level_start(self):
        """Prior median of the level, used as a starting value."""
        return -2.0 * np.log(np.log(2.0) / self.tau_level_rate)


# Log precisions above this value are treated as infinite smoothing: the
# likelihood is exactly flat beyond it, which keeps the Cholesky factors of the
# posterior precision well conditioned.
TAU_CAP = 30.0


def level_logprior(level, pri: Priors):
    """Log density of the level of the log precisions (up to a constant),
    including the Jacobian of ``sd = exp(-level / 2)``."""
    return -pri.tau_level_rate * np.exp(-0.5 * level) - 0.5 * level


def rw_logprior(tau, zeta, pri: Priors):
    """Log density of the mean-anchored random walk prior on the log precisions
    (up to a constant). The random walk part is intrinsic of rank ``len(tau) - 1``."""
    tau = np.asarray(tau, float)
    kt = tau.size
    return (-(kt - 1) * np.log(zeta)
            - np.sum(np.diff(tau) ** 2) / (2.0 * zeta**2)
            + level_logprior(tau.mean(), pri))


def pc_logprior_logzeta(logzeta, pri: Priors):
    """Log density of ``log(zeta)`` under the exponential (PC) prior on ``zeta``,
    including the Jacobian of the log transform."""
    zeta = np.exp(logzeta)
    return np.log(pri.zeta_rate) - pri.zeta_rate * zeta + logzeta


def simulate_tau(rng, kt, pri: Priors, zeta=None):
    """Draw ``tau`` from its prior (used by the validation tests)."""
    if zeta is None:
        zeta = rng.exponential(1.0 / pri.zeta_rate)
    level = -2.0 * np.log(rng.exponential(1.0 / pri.tau_level_rate))
    inc = rng.normal(size=kt - 1)
    walk = np.concatenate([[0.0], np.cumsum(inc)])
    return level + zeta * (walk - walk.mean()), zeta
