"""Lag matrices, B-spline bases, end constraints and sufficient statistics."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import null_space


def lag_matrix(x, p, start_at_zero=True):
    """Matrix of lagged values of a series.

    Row ``t`` holds ``x[t+p], x[t+p-1], ..., x[t]``, i.e. lags 0 to ``p`` of the
    observation at time ``t + p``; the first ``p`` time points, for which the
    full set of lags is unavailable, are dropped.

    Parameters
    ----------
    x : array_like
        The series to lag.
    p : int
        Maximum lag.
    start_at_zero : bool
        If ``True`` (default) the first column is lag 0 and there are ``p + 1``
        columns; otherwise columns start at lag 1 and there are ``p`` columns.
    """
    x = np.asarray(x, dtype=float).ravel()
    p = int(p)
    if p < 1:
        raise ValueError("p must be a positive integer")
    if p >= x.size:
        raise ValueError("p must be smaller than len(x)")
    n = x.size - p
    idx = np.arange(n)[:, None] + p - np.arange(p + 1)[None, :]
    m = x[idx]
    return m if start_at_zero else m[:, 1:]


def bspline_knots(lo, hi, nseg, deg):
    """Equally spaced knots with ``deg`` extra knots beyond each end."""
    dx = (hi - lo) / nseg
    return lo + np.arange(-deg, nseg + deg + 1) * dx


def bspline_basis(x, nseg=10, deg=3, domain=None, derivative=0):
    """Equally spaced B-spline basis (Eilers and Marx, 1996) evaluated at ``x``.

    Returns an array of shape ``(len(x), nseg + deg)``. With ``derivative=1``
    the first derivatives of the basis functions are returned instead.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    lo, hi = (float(x.min()), float(x.max())) if domain is None else map(float, domain)
    t = bspline_knots(lo, hi, nseg, deg)
    k = nseg + deg
    spl = BSpline(t, np.eye(k), deg, extrapolate=True)
    if derivative:
        spl = spl.derivative(derivative)
    return np.asarray(spl(x))


def difference_matrix(k, order=1):
    d = np.eye(k)
    for _ in range(order):
        d = np.diff(d, axis=0)
    return d


def cumsum_matrix(k):
    """``T`` such that ``T @ d`` is ``(0, d_1, d_1 + d_2, ...)`` for ``d`` of length ``k - 1``."""
    t = np.tril(np.ones((k, k - 1)), k=-1)
    return t


END_CONSTRAINTS = ("none", "zero", "zero_slope")


@dataclass
class Structure:
    """Everything about the lag-curve parameterisation that does not depend on data."""

    p: int
    K: int
    deg: int
    m: int                   # number of linear end constraints
    B: np.ndarray            # (p + 1, K) full basis
    C: np.ndarray            # (m, K) constraint rows
    Z: np.ndarray            # (K, K_r) orthonormal basis of the constraint null space
    B_r: np.ndarray          # (p + 1, K_r) reduced basis
    D_r: np.ndarray          # (K - 1, K_r) reduced first-difference operator
    T: np.ndarray            # (K, K - 1) cumulative-sum matrix for the non-centred form
    a_slope: np.ndarray | None   # (K - 1,) slope constraint in difference space (m == 2)
    logdet_const: float      # m <= 1: log|Q_r(tau)| = logdet_const + sum(tau)

    @property
    def K_r(self):
        return self.Z.shape[1]

    @property
    def K_tau(self):
        return self.K - 1

    def lam_matrix(self, tau):
        """Prior precision ``D_r' diag(exp(tau)) D_r`` of the reduced coefficients."""
        from .priors import TAU_CAP

        return (self.D_r.T * np.exp(np.minimum(tau, TAU_CAP))) @ self.D_r

    def logdet_prior(self, tau):
        """Generalised log-determinant of the prior precision as a function of tau."""
        from .priors import TAU_CAP

        tau = np.minimum(tau, TAU_CAP)
        if self.m <= 1:
            return self.logdet_const + float(np.sum(tau))
        return float(np.linalg.slogdet(self.lam_matrix(tau))[1])

    def coefficients_from_increments(self, d, level=None):
        """Full coefficient vector ``b`` from first differences ``d`` and, for
        unconstrained curves, the level ``b_1``; constrained curves determine the
        level from the value constraint."""
        b_raw = self.T @ d
        if self.m == 0:
            return level + b_raw
        return b_raw - self.C[0] @ b_raw


def build_structure(p, K, deg=3, end_constraint="zero"):
    if end_constraint not in END_CONSTRAINTS:
        raise ValueError(f"end_constraint must be one of {END_CONSTRAINTS}")
    if K < deg + 1:
        raise ValueError(f"K must be at least {deg + 1}")
    lags = np.arange(p + 1, dtype=float)
    B = bspline_basis(lags, nseg=K - deg, deg=deg, domain=(0.0, float(p)))
    rows = []
    if end_constraint in ("zero", "zero_slope"):
        rows.append(B[-1])
    if end_constraint == "zero_slope":
        rows.append(bspline_basis([float(p)], nseg=K - deg, deg=deg, domain=(0.0, float(p)), derivative=1)[0])
    C = np.array(rows).reshape(len(rows), K)
    m = len(rows)
    Z = np.eye(K) if m == 0 else null_space(C)
    D = difference_matrix(K)
    D_r = D @ Z
    T = cumsum_matrix(K)
    a_slope = (C[1] @ T) if m == 2 else None
    if m == 0:
        logdet_const = float(np.log(K))              # matrix-tree theorem: pdet = K prod(lambda)
    elif m == 1:
        logdet_const = float(np.linalg.slogdet(D_r.T @ D_r)[1])
    else:
        logdet_const = float("nan")
    return Structure(p=p, K=K, deg=deg, m=m, B=B, C=C, Z=Z, B_r=B @ Z, D_r=D_r, T=T,
                     a_slope=a_slope, logdet_const=logdet_const)


@dataclass
class SuffStats:
    """Sufficient statistics of the Gaussian likelihood as quadratics in the AR(1)
    coefficient ``phi``: ``S(phi) = S_a + phi S_b + phi^2 S_c`` and likewise for
    ``s = X'y`` and ``yy = y'y`` after prewhitening."""

    S_a: np.ndarray
    S_b: np.ndarray
    S_c: np.ndarray
    s_a: np.ndarray
    s_b: np.ndarray
    s_c: np.ndarray
    yy_a: float
    yy_b: float
    yy_c: float
    N: int

    @staticmethod
    def from_data(X, y):
        X = np.asarray(X, float)
        y = np.asarray(y, float)
        X0, X1, Xm = X[:-1], X[1:], X[1:-1]
        y0, y1, ym = y[:-1], y[1:], y[1:-1]
        return SuffStats(
            S_a=X.T @ X, S_b=-(X1.T @ X0 + X0.T @ X1), S_c=Xm.T @ Xm,
            s_a=X.T @ y, s_b=-(X1.T @ y0 + X0.T @ y1), s_c=Xm.T @ ym,
            yy_a=float(y @ y), yy_b=float(-2.0 * y1 @ y0), yy_c=float(ym @ ym), N=X.shape[0],
        )

    def at(self, phi=0.0):
        return (self.S_a + phi * self.S_b + phi**2 * self.S_c,
                self.s_a + phi * self.s_b + phi**2 * self.s_c,
                self.yy_a + phi * self.yy_b + phi**2 * self.yy_c)


def ar1_whitening_matrix(n, phi):
    """Explicit (n x n) matrix ``L`` with ``L' L`` the inverse correlation matrix
    of a stationary AR(1) process; used only in tests."""
    L = np.eye(n)
    L[0, 0] = np.sqrt(1.0 - phi**2)
    for t in range(1, n):
        L[t, t - 1] = -phi
    return L


@dataclass
class Standardizer:
    x_mean: float
    x_sd: float
    y_mean: float
    y_sd: float

    def x(self, x):
        return (np.asarray(x, float) - self.x_mean) / self.x_sd

    def y(self, y):
        return (np.asarray(y, float) - self.y_mean) / self.y_sd


@dataclass
class DesignData:
    """Data-dependent design for one fit (standardised scale)."""

    L: np.ndarray            # (N, p + 1) lag matrix of the standardised exposure
    W: np.ndarray            # (N, q) unpenalised columns (intercept, covariates)
    X: np.ndarray            # (N, K_r + q) [L B_r, W]
    y: np.ndarray            # (N,) response on the model scale
    offset: np.ndarray       # (N,)
    stats: SuffStats | None  # Gaussian sufficient statistics

    @property
    def N(self):
        return self.X.shape[0]

    @property
    def q(self):
        return self.W.shape[1]


def build_design(x, y, struct, covariates=None, intercept=True, offset=None,
                 standardize=True, family="gaussian", groups=None):
    """Assemble the design for one fit.

    ``groups`` (same length as ``x``) marks separate series, for example
    trading days or cities; lags are then formed within each group and the
    first ``p`` observations of every group are dropped. ``y``, ``covariates``
    and ``offset`` may be given either for all of ``x`` or only for the
    retained rows.
    """
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    p = struct.p
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError("x and y must be finite; missing values are not supported")
    if p < 1 or p >= x.size:
        raise ValueError("nlag must be between 1 and len(x) - 1")
    if groups is None:
        segments = [np.arange(x.size)]
    else:
        groups = np.asarray(groups).ravel()
        if groups.size != x.size:
            raise ValueError("groups must have the same length as x")
        _, first = np.unique(groups, return_index=True)
        segments = [np.where(groups == groups[i])[0] for i in np.sort(first)]
        if min(len(seg) for seg in segments) <= p:
            raise ValueError("every group must have more than nlag observations")
    keep = np.concatenate([seg[p:] for seg in segments])
    N = keep.size

    def rows(v, name):
        v = np.asarray(v, float)
        if v.shape[0] == x.size:
            return v[keep]
        if v.shape[0] == N:
            return v
        raise ValueError(f"{name} must have len(x) or len(x) - nlag * n_groups rows")

    y = rows(y, "y")
    if family == "gaussian" and standardize:
        y_mean = float(y.mean()) if intercept else 0.0
        y_sd = float(y.std(ddof=1)) if y.std(ddof=1) > 0 else 1.0
    else:
        y_mean, y_sd = 0.0, 1.0
    x_mean = float(x.mean()) if standardize else 0.0
    x_sd = float(x.std(ddof=1)) if (standardize and x.std(ddof=1) > 0) else 1.0
    std = Standardizer(x_mean=x_mean, x_sd=x_sd, y_mean=y_mean, y_sd=y_sd)

    xs = std.x(x)
    L = np.vstack([lag_matrix(xs[seg], p) for seg in segments])
    cols = []
    if intercept:
        cols.append(np.ones((N, 1)))
    if covariates is not None:
        Wc = np.asarray(covariates, float)
        Wc = rows(Wc.reshape(Wc.shape[0], -1), "covariates")
        cols.append(Wc)
    W = np.hstack(cols) if cols else np.zeros((N, 0))
    off = np.zeros(N) if offset is None else rows(np.asarray(offset, float).ravel(), "offset")
    X = np.hstack([L @ struct.B_r, W])
    y_model = std.y(y) if family == "gaussian" else y
    stats = SuffStats.from_data(X, y_model - off) if family == "gaussian" else None
    return DesignData(L=L, W=W, X=X, y=y_model, offset=off, stats=stats), std
