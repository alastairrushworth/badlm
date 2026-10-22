"""Walk-forward evaluation of a lag-curve forecaster and cost-aware trading.

:class:`FastFitter` compiles the collapsed posterior once for a given model
configuration and finds its mode for any set of sufficient statistics, so a
rolling re-fit costs an L-BFGS run from the previous mode rather than a fresh
compilation. :func:`walk_forward` refits at a chosen cadence and returns the
out-of-sample predictive mean and variance of every row; :func:`trade`
turns predictions into positions with proportional costs, either by trading a
fraction of the distance to the aim portfolio each period (Garleanu and
Pedersen, 2013) or with a no-trade band; :func:`performance` reports gross and
net Sharpe ratios, turnover and the break-even cost, and
:func:`deflated_sharpe` adjusts for the number of configurations tried
(Bailey and Lopez de Prado, 2014).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import jax
import jax.numpy as jnp
import numpy as np
from jax.flatten_util import ravel_pytree
from numpyro.infer import init_to_value
from numpyro.infer.util import initialize_model
from scipy.optimize import minimize
from scipy.special import ndtr, ndtri

from .collapsed import multi_precision
from .design import SuffStats, build_multi_design


# ---------------------------------------------------------------------------
# Fast mode finding with a reusable compiled objective
# ---------------------------------------------------------------------------
class FastFitter:
    """Posterior mode of the collapsed Gaussian model for arbitrary sufficient
    statistics, compiled once.

    Parameters
    ----------
    model : BADLM
        A fitted-or-not estimator whose ``components_`` are set (call
        ``prepare`` with design data first if they are not).
    stats : SuffStats
        Statistics used to trace the objective; later calls may pass any
        statistics of the same shapes.
    """

    def __init__(self, model, stats, q):
        from .numpyro_model import gaussian_model, ncp_init_values

        self.model = model
        self.q = q
        self.comps = model.components_
        start = ncp_init_values(self.comps, model.tau_model, model.priors, model.C_, model.ar1)
        args = (stats, self.comps, q, model.ar1, model.tau_model, model.priors, model.C_)
        info = initialize_model(jax.random.PRNGKey(0), gaussian_model, model_args=args,
                                init_strategy=init_to_value(values=start), dynamic_args=True)
        z0, self.unravel = ravel_pytree(info.param_info.z)
        self.z0 = np.asarray(z0, float)
        self.names = list(info.param_info.z.keys())
        comps, ar1, tm, pri, C = self.comps, model.ar1, model.tau_model, model.priors, model.C_

        def potential(zf, st):
            return info.potential_fn(st, comps, q, ar1, tm, pri, C)(self.unravel(zf))

        def constrain(zf, st):
            return info.postprocess_fn(st, comps, q, ar1, tm, pri, C)(self.unravel(zf))

        self._vg = jax.jit(jax.value_and_grad(potential))
        self._constrain = jax.jit(constrain)
        self._coef = jax.jit(self._coef_fn)
        self.last_mode = self.z0.copy()

    def _coef_fn(self, zf, st):
        site = self._constrain_impl(zf, st)
        taus = [site["tau" + ("" if e == 0 else f"_{e}")] for e in range(len(self.comps))]
        kappas = [site["mod_kappa" + ("" if e == 0 else f"_{e}")] if c.modifier is not None else None
                  for e, c in enumerate(self.comps)]
        Q, _ = multi_precision(self.comps, taus, kappas)
        sigma2 = site["sigma2"]
        phi = site["phi"] if self.model.ar1 else 0.0
        S = st.S_a + phi * st.S_b + phi**2 * st.S_c
        s = st.s_a + phi * st.s_b + phi**2 * st.s_c
        k = Q.shape[0]
        A = jnp.zeros((k + self.q, k + self.q)).at[:k, :k].set(Q) + S / sigma2
        A_inv = jnp.linalg.inv(A)
        mean = A_inv @ (s / sigma2)
        return mean, A_inv, sigma2, taus

    def _constrain_impl(self, zf, st):
        return self._constrain(zf, st)

    def fit(self, stats, warm=True, maxiter=500):
        """Mode of the posterior for ``stats``; returns the coefficient posterior
        mean, its covariance, the error variance and the mode's log precisions."""
        z_start = self.last_mode if warm else self.z0

        def f(z):
            v, g = self._vg(jnp.asarray(z), stats)
            v = float(v)
            g = np.asarray(g, float)
            if not np.isfinite(v) or not np.all(np.isfinite(g)):
                return 1e30, np.zeros_like(g)
            return v, g

        res = minimize(f, z_start, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-6})
        if not np.all(np.isfinite(res.x)):
            res = minimize(f, self.z0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
        self.last_mode = np.asarray(res.x, float)
        mean, cov, sigma2, taus = self._coef(jnp.asarray(self.last_mode), stats)
        return {"mean": np.asarray(mean), "cov": np.asarray(cov), "sigma2": float(sigma2),
                "taus": [np.asarray(t) for t in taus], "converged": bool(res.success), "nit": res.nit}


# ---------------------------------------------------------------------------
# Walk-forward predictions
# ---------------------------------------------------------------------------
@dataclass
class WalkForwardResult:
    time: np.ndarray            # time stamp of each predicted row (target time)
    row: np.ndarray             # index of the predicted row in the input arrays
    mean: np.ndarray            # predictive mean of the response (model scale)
    var: np.ndarray             # predictive variance (model scale)
    y: np.ndarray               # realised response
    refit_times: list = field(default_factory=list)
    curves: list = field(default_factory=list)      # posterior-mode lag curve(s) at each refit
    fit_info: list = field(default_factory=list)


def walk_forward(model, x, y, time, *, groups=None, covariates=None, weights=None, modifier=None,
                 start, step, min_train=None, window=None, verbose=False):
    """Expanding (or rolling) walk-forward fit.

    ``time`` is the time stamp of each row's response (the moment the
    predicted quantity is realised); a row is in the training set at a refit
    time ``T`` when ``time <= T`` and is predicted by the fit at ``T`` when
    ``T < time <= T + step``. ``start`` is the first refit time, ``step`` the
    refit cadence in the units of ``time``, ``window`` an optional rolling
    length. Inputs must already be on the scale the model should see
    (``standardize`` is switched off so no full-sample statistic leaks).

    The predictive mean is the regression mean ``X beta``. With ``ar1=True``
    the AR(1) model would also carry over ``phi`` times the previous residual,
    which this function does not add, so do not use AR(1) errors for
    prediction (with overlapping targets the regression mean alone is badly
    wrong).
    """
    model.standardize = False
    x_arr = np.asarray(x, float)
    if x_arr.ndim == 1:
        x_arr = x_arr[:, None]
    if modifier is not None and model.modifier_k is not None:
        from .design import LagComponent, build_modifier

        kind = modifier if isinstance(modifier, str) else "covariate"
        src = None if kind == "time" else (x_arr[:, 0] if kind == "exposure" else modifier)
        mod_basis = build_modifier(src, K=model.modifier_k, kind=kind)
    else:
        mod_basis = None
    from .design import LagComponent

    model.components_ = [LagComponent(st, mod_basis) for st in model.structs_]
    from .numpyro_model import tau_basis_matrix

    model.C_ = tau_basis_matrix(model.struct.K_tau, model.priors.tau_basis) if model.tau_model == "basis" else None
    data, _ = build_multi_design(x_arr, y, model.components_, covariates=covariates, intercept=model.intercept,
                                 standardize=False, groups=groups, weights=weights, modifier=modifier)
    t_row = np.asarray(time)[data.keep]
    X, yv = data.X, data.y
    w = data.weights
    fitter = None
    out_mean = np.full(yv.size, np.nan)
    out_var = np.full(yv.size, np.nan)
    refits, curves, infos = [], [], []
    T = start
    t_max = t_row.max()
    while T < t_max:
        train = t_row <= T
        if window is not None:
            train &= t_row > T - window
        test = (t_row > T) & (t_row <= T + step)
        if min_train is not None and train.sum() < min_train:
            T = T + step
            continue
        if train.sum() > X.shape[1] + 5 and test.any():
            stats = SuffStats.from_data(X[train], yv[train], weights=None if w is None else w[train])
            if fitter is None:
                fitter = FastFitter(model, stats, data.q)
            res = fitter.fit(stats)
            Xt = X[test]
            out_mean[test] = Xt @ res["mean"]
            out_var[test] = res["sigma2"] + np.einsum("ij,jk,ik->i", Xt, res["cov"], Xt)
            refits.append(T)
            curves.append(_mode_curves(model, res["mean"]))
            infos.append({"sigma2": res["sigma2"], "nit": res["nit"], "converged": res["converged"],
                          "n_train": int(train.sum())})
            if verbose:
                print(f"refit at {T}: n_train={train.sum()}, nit={res['nit']}, sigma2={res['sigma2']:.4f}")
        T = T + step
    have = np.isfinite(out_mean)
    return WalkForwardResult(time=t_row[have], row=data.keep[have], mean=out_mean[have], var=out_var[have],
                             y=yv[have], refit_times=refits, curves=curves, fit_info=infos)


def _mode_curves(model, theta):
    off, curves = 0, []
    for comp in model.components_:
        b = theta[off: off + comp.dim]
        off += comp.dim
        if comp.modifier is None:
            curves.append(comp.curve(b))
        else:
            curves.append(comp.curve(b, z=comp.modifier.grid(5)))
    return curves


# ---------------------------------------------------------------------------
# Trading and evaluation
# ---------------------------------------------------------------------------
def trade(mean, var, asset, time, *, vol=None, cost=0.0, rule="gp", speed=1.0, band=0.0,
          risk_aversion=1.0, clip=None):
    """Positions from predictive means and variances.

    Rows are (asset, time) pairs with the predicted response the next-period
    return in volatility-scaled units; ``vol`` (same rows) converts the aim
    position to notional so that each asset carries equal ex-ante risk. The aim
    is the Kelly/mean-variance position ``mean / (risk_aversion * var)``;
    ``rule='gp'`` moves a fraction ``speed`` of the way from the previous
    position to the aim each period, ``rule='band'`` trades to the aim only
    when the distance exceeds ``band`` (in aim units) and ``rule='aim'`` holds
    the aim. Returns a dict of arrays aligned with the inputs: ``position``
    (notional), ``trade`` (absolute notional traded).
    """
    mean, var, asset, time = map(np.asarray, (mean, var, asset, time))
    aim = mean / (risk_aversion * var)
    if clip is not None:
        aim = np.clip(aim, -clip, clip)
    if vol is not None:
        aim = aim / np.asarray(vol)
    pos = np.zeros_like(aim)
    trd = np.zeros_like(aim)
    for a in np.unique(asset):
        idx = np.where(asset == a)[0]
        idx = idx[np.argsort(time[idx], kind="stable")]
        prev = 0.0
        scale = np.mean(np.abs(aim[idx])) + 1e-12
        for i in idx:
            target = aim[i]
            if rule == "gp":
                new = prev + speed * (target - prev)
            elif rule == "band":
                new = target if abs(target - prev) > band * scale else prev
            else:
                new = target
            pos[i] = new
            trd[i] = abs(new - prev)
            prev = new
    return {"position": pos, "trade": trd, "aim": aim}


def performance(position, trade_size, realised, time, cost=0.0, periods_per_year=252, costs=(0.0,)):
    """Aggregate positions over assets by time and report performance.

    ``realised`` is the raw next-period return of each row (notional units),
    ``cost`` a proportional cost per unit of notional traded. Returns a dict
    with annualised gross and net Sharpe ratios for each cost level in
    ``costs``, annual turnover (notional traded / average gross notional),
    the break-even cost, and the daily portfolio return series.
    """
    position, trade_size, realised, time = map(np.asarray, (position, trade_size, realised, time))
    ts = np.unique(time)
    inv = np.searchsorted(ts, time)
    pnl = np.bincount(inv, weights=position * realised, minlength=ts.size)
    traded = np.bincount(inv, weights=trade_size, minlength=ts.size)
    gross_notional = np.bincount(inv, weights=np.abs(position), minlength=ts.size)
    out = {"time": ts, "pnl_gross": pnl, "traded": traded}
    sd = pnl.std(ddof=1)
    out["sharpe_gross"] = float(pnl.mean() / sd * np.sqrt(periods_per_year)) if sd > 0 else np.nan
    out["turnover_annual"] = float(traded.sum() / max(gross_notional.mean(), 1e-12) / ts.size * periods_per_year)
    out["breakeven_cost"] = float(pnl.sum() / max(traded.sum(), 1e-12))
    out["sharpe_net"] = {}
    for c in costs:
        net = pnl - c * traded
        s = net.std(ddof=1)
        out["sharpe_net"][c] = float(net.mean() / s * np.sqrt(periods_per_year)) if s > 0 else np.nan
    if cost:
        out["pnl_net"] = pnl - cost * traded
    return out


def probabilistic_sharpe(returns, benchmark_sr=0.0):
    """Probability that the true Sharpe ratio exceeds ``benchmark_sr`` (both per
    period), allowing for skewness and kurtosis of the returns."""
    r = np.asarray(returns, float)
    n = r.size
    sr = r.mean() / r.std(ddof=1)
    z = (r - r.mean()) / r.std(ddof=1)
    g3, g4 = np.mean(z**3), np.mean(z**4)
    denom = np.sqrt(max(1.0 - g3 * sr + (g4 - 1.0) / 4.0 * sr**2, 1e-12))
    return float(ndtr((sr - benchmark_sr) * np.sqrt(n - 1) / denom))


def deflated_sharpe(returns, n_trials, trial_sr_var=None):
    """Deflated Sharpe ratio of Bailey and Lopez de Prado (2014): the
    probabilistic Sharpe ratio against the expected maximum Sharpe ratio of
    ``n_trials`` unskilled strategies whose Sharpe ratios have variance
    ``trial_sr_var`` (per period). Returns ``(dsr, sr0)``."""
    r = np.asarray(returns, float)
    if trial_sr_var is None:
        trial_sr_var = (r.mean() / r.std(ddof=1)) ** 2 / 2.0
    euler = 0.5772156649
    n = max(int(n_trials), 2)
    sr0 = np.sqrt(trial_sr_var) * ((1 - euler) * ndtri(1 - 1.0 / n) + euler * ndtri(1 - 1.0 / (n * np.e)))
    return probabilistic_sharpe(r, sr0), float(sr0)
