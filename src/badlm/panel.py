"""Panel lag-curve models: every asset responds to a few driver series
through shared adaptive lag shapes and asset-specific loadings.

For asset ``i`` and drivers ``r = 1..R`` (a driver is a common series such as
a factor, or an asset-specific series such as the asset's own past),

    y_{i,t} = sum_r u_{ir} sum_l g_r(l) z_{r,i,t-l} + w_{i,t}' gamma + e_{i,t},

with ``g_r`` an adaptive lag curve (integrated out as in the single-series
model) and ``u_{ir} = 1 + s_r xi_{ir}``, ``xi ~ N(0, 1)``, ``s_r`` with a
penalised-complexity prior: the loadings shrink towards one, i.e. towards the
pooled model with one curve per driver, and ``s_r`` says how much assets
differ. With ``R`` drivers and ``N`` assets this is a rank-``R`` lag-curve
matrix ``beta_{ij}(l) = sum_r u_{ir} v_{jr} g_r(l)`` whose right factors
``v`` are fixed by the driver definitions.

The sufficient statistics of asset ``i`` enter as ``D_i S_i D_i`` and
``D_i s_i`` with ``D_i = diag(u_{i1} 1_K, ..., u_{iR} 1_K, 1_q)``, so the
collapsed likelihood is a cheap function of the loadings and NUTS or L-BFGS
can treat them as hyperparameters.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from jax.flatten_util import ravel_pytree
from numpyro.infer import init_to_value
from numpyro.infer.util import initialize_model
from scipy.optimize import minimize

from .collapsed import collapsed_loglik_Q, multi_precision
from .design import LagComponent, SuffStats, build_structure, lag_matrix
from .numpyro_model import _init_one, sample_kappa, sample_tau, suffix
from .priors import Priors


@dataclass
class PanelDesign:
    """Rows of a stacked panel with the driver lag design.

    ``X`` holds ``[L_1 B_1, ..., L_R B_R, W]`` for every retained row,
    ``asset`` the asset index of each row (0..N-1), ``time`` the row's time
    stamp and ``y`` the response.
    """

    X: np.ndarray
    y: np.ndarray
    asset: np.ndarray
    time: np.ndarray
    row: np.ndarray                 # index of the row in the input arrays
    comps: list
    q: int
    n_assets: int
    names: list = field(default_factory=list)

    @property
    def dims(self):
        return [c.dim for c in self.comps]

    def stats(self, mask):
        """Per-asset sufficient statistics over the rows selected by ``mask``:
        arrays ``S`` (N, d, d), ``s`` (N, d), ``yy`` (N,), ``n`` (N,)."""
        d = self.X.shape[1]
        S = np.zeros((self.n_assets, d, d)); s = np.zeros((self.n_assets, d))
        yy = np.zeros(self.n_assets); n = np.zeros(self.n_assets)
        for i in range(self.n_assets):
            m = mask & (self.asset == i)
            if m.any():
                Xi, yi = self.X[m], self.y[m]
                S[i] = Xi.T @ Xi; s[i] = Xi.T @ yi; yy[i] = yi @ yi; n[i] = m.sum()
        return {"S": S, "s": s, "yy": yy, "n": n}


def build_panel_design(drivers, y, time, asset, nlag, k=None, end_constraint="zero", intercept=True,
                       lag_scale="linear", standardize_drivers=True, modifier=None, modifier_k=5):
    """Assemble a :class:`PanelDesign`.

    ``drivers`` is a list of ``R`` arrays, each aligned with ``y`` (one value
    per row of the stacked panel; a common factor is simply repeated on every
    asset's rows). Lags are formed within each asset's rows, which must be
    contiguous and time-ordered. ``nlag`` and ``k`` are scalars or one per
    driver. ``modifier`` (aligned with ``y``) makes every driver's curve a
    surface over lag and that state, as in :class:`badlm.BADLM`.
    """
    from .design import build_modifier

    mod = None if modifier is None else build_modifier(np.asarray(modifier, float), K=modifier_k)
    y = np.asarray(y, float).ravel()
    asset = np.asarray(asset)
    names, inv = np.unique(asset, return_inverse=True)
    R = len(drivers)
    nlags = [int(nlag)] * R if np.ndim(nlag) == 0 else [int(v) for v in nlag]
    ks = [None] * R if k is None else ([int(k)] * R if np.ndim(k) == 0 else [int(v) for v in k])
    comps = []
    for p, kk in zip(nlags, ks):
        kk = int(np.ceil(2.0 * (p + 1) / 3.0)) if kk is None else kk
        comps.append(LagComponent(build_structure(p, kk, end_constraint=end_constraint, lag_scale=lag_scale), mod))
    p_max = max(nlags)
    Z = [np.asarray(d, float).ravel() for d in drivers]
    if standardize_drivers:
        Z = [(z - z.mean()) / (z.std(ddof=1) if z.std(ddof=1) > 0 else 1.0) for z in Z]
    blocks, keep_rows, keep_asset = [[] for _ in range(R)], [], []
    for i in range(names.size):
        rows = np.where(inv == i)[0]
        if rows.size <= p_max:
            continue
        z_rows = None if mod is None else np.asarray(modifier, float)[rows[p_max:]]
        for r in range(R):
            L = lag_matrix(Z[r][rows], nlags[r])[p_max - nlags[r]:]
            blocks[r].append(comps[r].design(L, z_rows))
        keep_rows.append(rows[p_max:])
        keep_asset.append(np.full(rows.size - p_max, i))
    keep_rows = np.concatenate(keep_rows)
    cols = [np.vstack(b) for b in blocks]
    if intercept:
        cols.append(np.ones((keep_rows.size, 1)))
    X = np.hstack(cols)
    return PanelDesign(X=X, y=y[keep_rows], asset=np.concatenate(keep_asset), time=np.asarray(time)[keep_rows],
                       row=keep_rows, comps=comps, q=int(intercept), n_assets=names.size, names=list(names))


def loading_matrix(u, dims, q):
    """Per-asset scaling vectors ``d_i`` (N, sum(dims) + q) from loadings ``u`` (N, R)."""
    parts = [jnp.repeat(u[:, r: r + 1], d, axis=1) for r, d in enumerate(dims)]
    if q:
        parts.append(jnp.ones((u.shape[0], q)))
    return jnp.concatenate(parts, axis=1)


def panel_stats(stats, d):
    """Total sufficient statistics for scaling vectors ``d`` (N, dim)."""
    S = jnp.einsum("ia,iab,ib->ab", d, stats["S"], d)
    s = jnp.einsum("ia,ia->a", d, stats["s"])
    yy = jnp.sum(stats["yy"])
    N = jnp.sum(stats["n"])
    zero = jnp.zeros_like(S)
    return SuffStats(S, zero, zero, s, jnp.zeros_like(s), jnp.zeros_like(s), yy, 0.0, 0.0, N)


def panel_model(stats, comps, q, n_assets, tau_model, pri: Priors, C=None, loading_scale=0.5, pool=False):
    """Collapsed panel model: curves integrated out, loadings and precisions
    sampled. ``pool=True`` fixes every loading at one (the pooled model)."""
    R = len(comps)
    taus = [sample_tau(c.n_tau, tau_model, pri, C, suffix(r)) for r, c in enumerate(comps)]
    kappas = [sample_kappa(pri, suffix(r)) if c.modifier is not None else None for r, c in enumerate(comps)]
    if pool:
        u = jnp.ones((n_assets, R))
    else:
        s_u = numpyro.sample("loading_sd", dist.Exponential(-np.log(0.1) / loading_scale).expand([R]))
        xi = numpyro.sample("loading_z", dist.Normal(0.0, 1.0).expand([n_assets, R]))
        u = numpyro.deterministic("loading", 1.0 + s_u[None, :] * xi)
    sigma2 = numpyro.sample("sigma2", dist.InverseGamma(pri.sigma2_a, pri.sigma2_b))
    d = loading_matrix(u, [c.dim for c in comps], q)
    tot = panel_stats(stats, d)
    Q, logdet = multi_precision(comps, taus, kappas)
    numpyro.factor("loglik", collapsed_loglik_Q(Q, logdet, sigma2, 0.0, tot, q, False) / pri.temper)


class PanelFitter:
    """Posterior mode of the panel model for any per-asset statistics of the
    same shapes, compiled once (see :class:`badlm.backtest.FastFitter`)."""

    def __init__(self, design: PanelDesign, stats, *, tau_model="rw", priors=None, loading_scale=0.5, pool=False):
        self.design = design
        self.pri = priors or Priors()
        self.pool = pool
        comps, q, N = design.comps, design.q, design.n_assets
        start = {"sigma2": 0.5}
        for r, c in enumerate(comps):
            start.update(_init_one(c.n_tau, tau_model, self.pri, None, suffix(r)))
            if c.modifier is not None:
                start["mod_sd" + suffix(r)] = 0.3
        if not pool:
            start.update({"loading_sd": np.full(len(comps), 0.1), "loading_z": np.zeros((N, len(comps)))})
        args = (stats, comps, q, N, tau_model, self.pri, None, loading_scale, pool)
        info = initialize_model(jax.random.PRNGKey(0), panel_model, model_args=args,
                                init_strategy=init_to_value(values=start), dynamic_args=True)
        z0, self.unravel = ravel_pytree(info.param_info.z)
        self.z0 = np.asarray(z0, float)
        rest = (comps, q, N, tau_model, self.pri, None, loading_scale, pool)

        def potential(zf, st):
            return info.potential_fn(st, *rest)(self.unravel(zf))

        def constrain(zf, st):
            return info.postprocess_fn(st, *rest)(self.unravel(zf))

        self._vg = jax.jit(jax.value_and_grad(potential))
        self._constrain = jax.jit(constrain)
        self.comps, self.q, self.N = comps, q, N
        self.last_mode = self.z0.copy()

    def fit(self, stats, warm=True, maxiter=500):
        z_start = self.last_mode if warm else self.z0

        def f(z):
            v, g = self._vg(jnp.asarray(z), stats)
            v, g = float(v), np.asarray(g, float)
            return (v, g) if np.isfinite(v) and np.all(np.isfinite(g)) else (1e30, np.zeros_like(g))

        res = minimize(f, z_start, jac=True, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-6})
        self.last_mode = np.asarray(res.x, float)
        site = {k: np.asarray(v) for k, v in self._constrain(jnp.asarray(self.last_mode), stats).items()}
        R = len(self.comps)
        taus = [site["tau" + suffix(r)] for r in range(R)]
        kappas = [jnp.asarray(site["mod_kappa" + suffix(r)]) if c.modifier is not None else None
                  for r, c in enumerate(self.comps)]
        u = np.ones((self.N, R)) if self.pool else site["loading"]
        d = np.asarray(loading_matrix(jnp.asarray(u), [c.dim for c in self.comps], self.q))
        tot = panel_stats(stats, jnp.asarray(d))
        Q, _ = multi_precision(self.comps, [jnp.asarray(t) for t in taus], kappas)
        k = Q.shape[0]
        sigma2 = float(site["sigma2"])
        A = np.zeros((k + self.q, k + self.q)); A[:k, :k] = np.asarray(Q)
        A = A + np.asarray(tot.S_a) / sigma2
        cov = np.linalg.inv(A)
        mean = cov @ (np.asarray(tot.s_a) / sigma2)
        return {"mean": mean, "cov": cov, "sigma2": sigma2, "loading": u, "d": d, "taus": taus,
                "loading_sd": None if self.pool else site["loading_sd"], "converged": bool(res.success), "nit": res.nit}

    def curves(self, res):
        """Posterior-mode lag curve of each driver for the average asset (for
        a modifier model, the curves at five values of the state)."""
        out, off = [], 0
        for c in self.comps:
            b = res["mean"][off: off + c.dim]
            out.append(c.curve(b) if c.modifier is None else c.curve(b, z=c.modifier.grid(5)))
            off += c.dim
        return out


def panel_walk_forward(design: PanelDesign, *, start, step, tau_model="rw", priors=None, loading_scale=0.5,
                       pool=False, min_train=None, verbose=False):
    """Expanding walk-forward fit of the panel model (see
    :func:`badlm.backtest.walk_forward` for the conventions). Returns a dict
    with out-of-sample ``mean``, ``var``, ``y``, ``time``, ``asset``, ``row``
    and the fitted loadings and curves at each refit."""
    t = design.time
    out_mean = np.full(design.y.size, np.nan); out_var = np.full(design.y.size, np.nan)
    fitter, refits = None, []
    T, t_max = start, t.max()
    while T < t_max:
        train = t <= T
        test = (t > T) & (t <= T + step)
        if (min_train is None or train.sum() >= min_train) and train.sum() > design.X.shape[1] + 5 and test.any():
            st = design.stats(train)
            if fitter is None:
                fitter = PanelFitter(design, st, tau_model=tau_model, priors=priors, loading_scale=loading_scale, pool=pool)
            res = fitter.fit(st)
            Xt = design.X[test] * res["d"][design.asset[test]]
            out_mean[test] = Xt @ res["mean"]
            out_var[test] = res["sigma2"] + np.einsum("ij,jk,ik->i", Xt, res["cov"], Xt)
            refits.append({"time": T, "loading": res["loading"], "loading_sd": res["loading_sd"],
                           "curves": fitter.curves(res), "sigma2": res["sigma2"], "nit": res["nit"]})
            if verbose:
                print(f"refit at {T}: n_train={train.sum()}, nit={res['nit']}, loading sd={res['loading_sd']}")
        T = T + step
    have = np.isfinite(out_mean)
    return {"mean": out_mean[have], "var": out_var[have], "y": design.y[have], "time": t[have],
            "asset": design.asset[have], "row": design.row[have], "refits": refits}
