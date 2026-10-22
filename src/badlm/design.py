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


LAG_SCALES = ("linear", "log")


def build_structure(p, K, deg=3, end_constraint="zero", lag_scale="linear"):
    """Lag-curve parameterisation. ``lag_scale='log'`` places the equally spaced
    B-spline knots on ``log(1 + lag)``, so that the curve can change quickly at
    short lags and only slowly at long ones."""
    if end_constraint not in END_CONSTRAINTS:
        raise ValueError(f"end_constraint must be one of {END_CONSTRAINTS}")
    if lag_scale not in LAG_SCALES:
        raise ValueError(f"lag_scale must be one of {LAG_SCALES}")
    if K < deg + 1:
        raise ValueError(f"K must be at least {deg + 1}")
    lags = np.arange(p + 1, dtype=float)
    u = np.log1p(lags) if lag_scale == "log" else lags
    dom = (0.0, float(u[-1]))
    B = bspline_basis(u, nseg=K - deg, deg=deg, domain=dom)
    rows = []
    if end_constraint in ("zero", "zero_slope"):
        rows.append(B[-1])
    if end_constraint == "zero_slope":
        rows.append(bspline_basis([u[-1]], nseg=K - deg, deg=deg, domain=dom, derivative=1)[0])
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
    def from_data(X, y, weights=None):
        """Sufficient statistics; with ``weights`` each row of ``X`` and ``y`` is
        multiplied by the square root of its weight (inverse-variance weighting)."""
        X = np.asarray(X, float)
        y = np.asarray(y, float)
        if weights is not None:
            sw = np.sqrt(np.asarray(weights, float).ravel())
            X = X * sw[:, None]
            y = y * sw
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
    x_mean: float | np.ndarray
    x_sd: float | np.ndarray
    y_mean: float
    y_sd: float

    def x(self, x):
        x = np.asarray(x, float)
        if x.ndim == 2 and np.ndim(self.x_mean) == 1:
            return (x - np.asarray(self.x_mean)[None, :]) / np.asarray(self.x_sd)[None, :]
        return (x - self.x_mean) / self.x_sd

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
    Ls: list | None = None   # lag matrices of every exposure (multi-exposure fits)
    z: np.ndarray | None = None        # modifier values per retained row
    weights: np.ndarray | None = None  # observation weights per retained row
    keep: np.ndarray | None = None     # indices of the retained rows in the input

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


# ----------------------------------------------------------------------------
# Modifier bases and lag components (multi-exposure and tensor-product curves)
# ----------------------------------------------------------------------------
@dataclass
class ModifierBasis:
    """B-spline basis over a modifier variable ``z`` (time, volatility, the
    exposure level) that lets the lag curve vary smoothly with ``z``.

    The roughness penalty along ``z`` is the first-difference penalty
    ``P = D'D`` whose null space is 'constant in z', i.e. the static curve, so
    shrinking its scale recovers the ordinary model.
    """

    K: int
    deg: int
    domain: tuple
    kind: str                    # 'covariate', 'time' or 'exposure'
    P: np.ndarray                # (K, K) first-difference penalty
    mu: np.ndarray               # eigenvalues of P (ascending; mu[0] == 0)
    V: np.ndarray                # eigenvectors of P

    def basis(self, z):
        z = np.clip(np.asarray(z, float).ravel(), self.domain[0], self.domain[1])
        return bspline_basis(z, nseg=self.K - self.deg, deg=self.deg, domain=self.domain)

    def grid(self, n=5):
        return np.linspace(self.domain[0], self.domain[1], n)


def build_modifier(z, K=6, deg=3, kind="covariate", domain=None, trim=0.01):
    """Modifier basis with the domain taken from the central ``1 - 2 trim`` of
    the training values (``kind='covariate'`` or ``'exposure'``) or ``[0, 1]``
    (``kind='time'``). With ``kind='exposure'`` the modifier of lag ``j`` is the
    lagged exposure itself, so the surface ``beta(j, x)`` is a non-linear
    exposure-lag response (the cross-basis of Gasparrini et al., 2010)."""
    if K < deg + 1:
        raise ValueError(f"modifier K must be at least {deg + 1}")
    if domain is None:
        if kind == "time":
            domain = (0.0, 1.0)
        else:
            z = np.asarray(z, float).ravel()
            lo, hi = np.quantile(z, [trim, 1.0 - trim])
            if hi <= lo:
                lo, hi = float(z.min()), float(z.max() + 1e-12)
            domain = (float(lo), float(hi))
    D = difference_matrix(K)
    P = D.T @ D
    mu, V = np.linalg.eigh(P)
    mu = np.maximum(mu, 0.0)
    mu[0] = 0.0
    return ModifierBasis(K=K, deg=deg, domain=tuple(domain), kind=kind, P=P, mu=mu, V=V)


@dataclass
class LagComponent:
    """One exposure's lag curve: a :class:`Structure` and, optionally, a
    :class:`ModifierBasis` that makes the curve a surface over lag and ``z``.

    Coefficients are ordered ``(k, l) -> k * K_z + l``: lag basis function ``k``,
    modifier basis function ``l``. The prior precision is the Kronecker sum
    ``Q_lag(tau) (x) I + kappa I (x) P_z`` whose log-determinant follows from the
    eigenvalues of the two factors.
    """

    struct: Structure
    modifier: ModifierBasis | None = None

    @property
    def K_z(self):
        return 1 if self.modifier is None else self.modifier.K

    @property
    def dim(self):
        return self.struct.K_r * self.K_z

    @property
    def n_tau(self):
        return self.struct.K_tau

    def design(self, L, z=None):
        """Design block for lag matrix ``L`` (standardised exposure) and modifier
        values ``z`` (one per row of ``L``). For an exposure modifier the block
        is ``sum_j x_{t-j} B(j) (x) C(x_{t-j})``."""
        if self.modifier is None:
            return L @ self.struct.B_r
        if self.modifier.kind == "exposure":
            B_r = self.struct.B_r
            out = np.zeros((L.shape[0], self.dim))
            for j in range(L.shape[1]):
                Cj = self.modifier.basis(L[:, j])                       # (N, K_z)
                out += ((L[:, j, None] * B_r[j][None, :])[:, :, None] * Cj[:, None, :]).reshape(L.shape[0], -1)
            return out
        LB = L @ self.struct.B_r
        Cz = self.modifier.basis(z)
        return (LB[:, :, None] * Cz[:, None, :]).reshape(L.shape[0], -1)

    def curve(self, b, z=None, weights=None):
        """Lag curve(s) from coefficients ``b`` (``..., dim``). Without a modifier
        the result has shape ``(..., p + 1)``. With one, ``z`` (scalar or array)
        gives the curve at those modifier values, shape ``(..., len(z), p + 1)``;
        alternatively ``weights`` (``K_z``,) averages the modifier basis, which
        gives the curve at the average training state."""
        b = np.asarray(b, float)
        B_r = self.struct.B_r
        if self.modifier is None:
            return b @ B_r.T
        bm = b.reshape(b.shape[:-1] + (self.struct.K_r, self.K_z))
        if weights is not None:
            return (bm @ np.asarray(weights, float)) @ B_r.T
        Cz = self.modifier.basis(np.atleast_1d(z))                # (nz, K_z)
        surf = np.einsum("...kl,zl->...zk", bm, Cz)
        return surf @ B_r.T


def build_multi_design(X_exp, y, components, covariates=None, intercept=True, offset=None,
                       standardize=True, family="gaussian", groups=None, weights=None, modifier=None):
    """Assemble the design for several exposures, optional modifier and weights.

    ``X_exp`` has shape ``(n, E)``; component ``e`` receives column ``e`` lagged
    to ``components[e].struct.p``. All lag matrices are aligned on the rows
    that have every lag available for the longest ``p``. ``modifier`` is a
    series of length ``n`` (or the string ``'time'``). ``weights`` (length ``n``
    or number of retained rows) multiply the squared residuals; they are the
    inverse of an ex-ante variance when the response is heteroskedastic.
    """
    X_exp = np.asarray(X_exp, float)
    if X_exp.ndim == 1:
        X_exp = X_exp[:, None]
    n, E = X_exp.shape
    if E != len(components):
        raise ValueError("number of exposure columns must match the number of components")
    y = np.asarray(y, float).ravel()
    if np.any(~np.isfinite(X_exp)) or np.any(~np.isfinite(y)):
        raise ValueError("x and y must be finite; missing values are not supported")
    p_max = max(c.struct.p for c in components)
    if p_max < 1 or p_max >= n:
        raise ValueError("nlag must be between 1 and len(x) - 1")
    if groups is None:
        segments = [np.arange(n)]
    else:
        groups = np.asarray(groups).ravel()
        if groups.size != n:
            raise ValueError("groups must have the same length as x")
        _, first = np.unique(groups, return_index=True)
        segments = [np.where(groups == groups[i])[0] for i in np.sort(first)]
        if min(len(seg) for seg in segments) <= p_max:
            raise ValueError("every group must have more than nlag observations")
    keep = np.concatenate([seg[p_max:] for seg in segments])
    N = keep.size

    def rows(v, name):
        v = np.asarray(v, float)
        if v.shape[0] == n:
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
    if standardize:
        x_mean = X_exp.mean(axis=0)
        x_sd = X_exp.std(axis=0, ddof=1)
        x_sd = np.where(x_sd > 0, x_sd, 1.0)
    else:
        x_mean, x_sd = np.zeros(E), np.ones(E)
    std = Standardizer(x_mean=x_mean, x_sd=x_sd, y_mean=y_mean, y_sd=y_sd)
    xs = std.x(X_exp)

    z = None
    if any(c.modifier is not None and c.modifier.kind != "exposure" for c in components):
        if modifier is None:
            raise ValueError("this model has a modifier; pass modifier= to fit")
        if isinstance(modifier, str) and modifier == "time":
            # position of each retained row in the input series, scaled to [0, 1];
            # for stacked panels pass calendar time as a numeric modifier instead
            z = (keep - keep.min()) / max(float(keep.max() - keep.min()), 1.0)
        else:
            zm = np.asarray(modifier, float)
            if zm.ndim == 2 and zm.shape[1] == E and E > 1:
                z = np.stack([rows(zm[:, e], "modifier") for e in range(E)], axis=1)   # one state per exposure
            else:
                z = rows(zm.ravel(), "modifier")

    blocks, Ls = [], []
    for e, comp in enumerate(components):
        p = comp.struct.p
        L = np.vstack([lag_matrix(xs[seg, e], p)[p_max - p:] for seg in segments])
        Ls.append(L)
        blocks.append(comp.design(L, z[:, e] if (z is not None and z.ndim == 2) else z))
    cols = []
    if intercept:
        cols.append(np.ones((N, 1)))
    if covariates is not None:
        Wc = np.asarray(covariates, float)
        Wc = rows(Wc.reshape(Wc.shape[0], -1), "covariates")
        cols.append(Wc)
    W = np.hstack(cols) if cols else np.zeros((N, 0))
    off = np.zeros(N) if offset is None else rows(np.asarray(offset, float).ravel(), "offset")
    X = np.hstack(blocks + [W])
    y_model = std.y(y) if family == "gaussian" else y
    w = None if weights is None else rows(np.asarray(weights, float).ravel(), "weights")
    if w is not None and np.any(w <= 0):
        raise ValueError("weights must be positive")
    stats = SuffStats.from_data(X, y_model - off, weights=w) if family == "gaussian" else None
    data = DesignData(L=Ls[0], W=W, X=X, y=y_model, offset=off, stats=stats, Ls=Ls, z=z,
                      weights=w, keep=keep)
    return data, std
