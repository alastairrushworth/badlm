"""Collapsed Gaussian likelihood: the spline coefficients (and any unpenalised
coefficients, which carry a flat prior) are integrated out analytically.

With ``theta = (b_r, gamma)``, prior precision ``Q = blockdiag(Q_r(tau), 0)`` and
whitened sufficient statistics ``S = X'X``, ``s = X'y``, ``yy = y'y``,

    log p(y | tau, sigma2, phi) = -N/2 log(2 pi sigma2) + q/2 log(2 pi)
        + 1/2 log|Q_r|  - 1/2 log|A| + 1/2 m' A^{-1} m - yy / (2 sigma2)
        [+ 1/2 log(1 - phi^2) for AR(1) errors]

where ``A = Q + S / sigma2`` and ``m = s / sigma2``. ``log|Q_r|`` is the
generalised determinant when ``Q_r`` is singular (no end constraint), which for
the chain penalty has the closed form ``log K + sum(tau)``.
"""
from __future__ import annotations

import jax.numpy as jnp
from jax.scipy.linalg import cho_solve

from .priors import TAU_CAP


def suff_at(stats, phi):
    S = stats.S_a + phi * stats.S_b + phi**2 * stats.S_c
    s = stats.s_a + phi * stats.s_b + phi**2 * stats.s_c
    yy = stats.yy_a + phi * stats.yy_b + phi**2 * stats.yy_c
    return S, s, yy


def logdet_prior(tau, struct):
    """Generalised log-determinant of the adaptive prior precision: closed form
    for at most one end constraint, Cholesky otherwise."""
    if struct.m <= 1:
        return struct.logdet_const + jnp.sum(tau)
    lam = jnp.exp(tau)
    Q_r = (jnp.asarray(struct.D_r).T * lam) @ jnp.asarray(struct.D_r)
    R = jnp.linalg.cholesky(Q_r)
    return 2.0 * jnp.sum(jnp.log(jnp.diag(R)))


def _embed(Q_r, q):
    if q == 0:
        return Q_r
    k = Q_r.shape[0]
    return jnp.zeros((k + q, k + q)).at[:k, :k].set(Q_r)


def collapsed_loglik_Q(Q_r, logdet_Q, sigma2, phi, stats, q, ar1=False):
    """Log marginal likelihood for an arbitrary prior precision ``Q_r`` of the
    penalised coefficients with (generalised) log-determinant ``logdet_Q``."""
    S, s, yy = suff_at(stats, phi)
    N = stats.N
    A = _embed(Q_r, q) + S / sigma2
    R = jnp.linalg.cholesky(A)
    logdet_A = 2.0 * jnp.sum(jnp.log(jnp.diag(R)))
    m = s / sigma2
    mu = cho_solve((R, True), m)
    ll = (-0.5 * N * jnp.log(2.0 * jnp.pi * sigma2)
          + 0.5 * q * jnp.log(2.0 * jnp.pi)
          + 0.5 * logdet_Q
          - 0.5 * logdet_A
          + 0.5 * jnp.dot(m, mu)
          - 0.5 * yy / sigma2)
    if ar1:
        ll = ll + 0.5 * jnp.log1p(-phi**2)
    return ll


def conditional_coefficients_Q(Q_r, sigma2, phi, stats, q, z):
    """Draw ``theta | Q_r, sigma2, phi, y`` given standard normal ``z`` and return
    the effective degrees of freedom ``tr[(S + sigma2 Q)^{-1} S]``."""
    S, s, _ = suff_at(stats, phi)
    A = _embed(Q_r, q) + S / sigma2
    R = jnp.linalg.cholesky(A)
    mu = cho_solve((R, True), s / sigma2)
    theta = mu + jnp.linalg.solve(R.T, z)
    ed = jnp.trace(cho_solve((R, True), S)) / sigma2
    return theta, ed


def adaptive_precision(tau, struct):
    """Reduced prior precision ``D_r' diag(exp(tau)) D_r`` of the adaptive model
    and its generalised log-determinant; log precisions are capped."""
    tau = jnp.minimum(tau, TAU_CAP)
    D_r = jnp.asarray(struct.D_r)
    return (D_r.T * jnp.exp(tau)) @ D_r, logdet_prior(tau, struct)


def collapsed_loglik(tau, sigma2, phi, stats, struct, q, ar1=False):
    """Log marginal likelihood of the adaptive model given its hyperparameters.
    Log precisions above ``TAU_CAP`` act as infinite smoothing: the likelihood
    is exactly flat in them, so the prior alone governs that region."""
    Q_r, logdet_Q = adaptive_precision(tau, struct)
    return collapsed_loglik_Q(Q_r, logdet_Q, sigma2, phi, stats, q, ar1)


def conditional_coefficients(tau, sigma2, phi, stats, struct, q, z):
    """Draw ``theta | tau, sigma2, phi, y`` for the adaptive model."""
    Q_r, _ = adaptive_precision(tau, struct)
    return conditional_coefficients_Q(Q_r, sigma2, phi, stats, q, z)
