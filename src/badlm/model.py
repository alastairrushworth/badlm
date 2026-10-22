"""User-facing estimator."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .design import (DesignData, LagComponent, Standardizer, Structure, build_modifier,
                     build_multi_design, build_structure, lag_matrix)
from .diagnostics import summarise_draws
from .priors import Priors

FAMILIES = ("gaussian", "poisson", "negbin")
SAMPLERS = ("nuts", "gibbs", "laplace")


@dataclass
class LagCurveSummary:
    lag: np.ndarray
    median: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    mean: np.ndarray
    probs: tuple

    def to_dataframe(self):
        import pandas as pd

        return pd.DataFrame({"lag": self.lag, "mean": self.mean, "median": self.median,
                             f"q{self.probs[0]:g}": self.lower, f"q{self.probs[1]:g}": self.upper})

    def __repr__(self):
        head = f"{'lag':>4} {'median':>10} {'lower':>10} {'upper':>10}\n"
        rows = [f"{int(l):>4} {m:>10.4g} {lo:>10.4g} {up:>10.4g}"
                for l, m, lo, up in zip(self.lag, self.median, self.lower, self.upper)]
        return head + "\n".join(rows)


@dataclass
class BADLM:
    """Bayesian adaptive distributed lag model.

    Parameters
    ----------
    nlag : int or sequence of int
        Maximum lag ``p``. With adaptive smoothing it can be chosen generously.
        A sequence gives one lag curve per exposure (``x`` then has one column
        per exposure); each curve has its own adaptive smoothing.
    k : int, sequence of int or None
        Number of cubic B-spline functions for the lag curve(s). Defaults to
        ``ceil(2 (p + 1) / 3)`` as in the paper's simulations.
    tau_model : {'rw', 'basis', 'monotone', 'horseshoe', 'constant'}
        Prior for the log precisions ``tau`` of successive differences of the
        spline coefficients. ``'rw'`` is the random walk of the paper with a
        penalised-complexity prior on its scale; ``'basis'`` represents ``tau``
        with a few B-splines; ``'monotone'`` forces smoothing to increase with
        lag; ``'horseshoe'`` uses independent local shrinkage; ``'constant'`` is
        the non-adaptive model (M1 in the paper).
    end_constraint : {'zero', 'zero_slope', 'none'}
        Exact linear constraint on the curve at lag ``p``, imposed by projecting
        onto the constraint null space.
    lag_scale : {'linear', 'log'}
        Knot placement of the lag basis: equally spaced in lag, or in
        ``log(1 + lag)`` so that short lags get more resolution than long ones.
    family : {'gaussian', 'poisson', 'negbin'}
    ar1 : bool
        AR(1) errors (Gaussian family, NUTS and Laplace only).
    intercept : bool
    standardize : bool
        Standardise ``x`` (and ``y`` for the Gaussian family) internally; all
        reported quantities are on the original scale.
    sampler : {'nuts', 'gibbs', 'laplace'}
        ``'nuts'`` integrates the coefficients out and runs NUTS on the
        hyperparameters; ``'gibbs'`` is a Gibbs/slice sampler (Gaussian, 'rw'
        only); ``'laplace'`` is a fast Gaussian approximation at the mode.
    priors : Priors
    parameterization : {'centred', 'noncentred'}
        For count families only: sample the spline coefficients directly
        (``'centred'``, good when the data are informative) or as scaled
        standard-normal increments (``'noncentred'``).
    modifier_k : int or None
        If given, every lag curve becomes a surface over lag and a modifier
        variable ``z`` (passed to ``fit`` as ``modifier=``) represented by
        ``modifier_k`` B-splines in ``z``, with the adaptive penalty along lag
        and a penalised-complexity prior on roughness along ``z`` that shrinks
        towards a curve that does not depend on ``z``. The modifier can be a
        state variable (volatility, volume), the exposure level, or
        ``'time'`` for a smoothly drifting curve. Gaussian family only.
    """

    nlag: int | tuple | list
    k: int | tuple | list | None = None
    tau_model: str = "rw"
    end_constraint: str = "zero"
    lag_scale: str = "linear"
    family: str = "gaussian"
    ar1: bool = False
    intercept: bool = True
    standardize: bool = True
    sampler: str = "nuts"
    priors: Priors = field(default_factory=Priors)
    parameterization: str = "centred"
    modifier_k: int | None = None

    def __post_init__(self):
        nlags = [int(self.nlag)] if np.ndim(self.nlag) == 0 else [int(v) for v in self.nlag]
        if self.k is None:
            ks = [int(np.ceil(2.0 * (p + 1) / 3.0)) for p in nlags]
        elif np.ndim(self.k) == 0:
            ks = [int(self.k)] * len(nlags)
        else:
            ks = [int(v) for v in self.k]
            if len(ks) != len(nlags):
                raise ValueError("k must have one entry per exposure")
        self.nlags_ = nlags
        self.ks_ = ks
        self.n_exposures = len(nlags)
        if self.n_exposures == 1:
            self.nlag, self.k = nlags[0], ks[0]
        multi = self.n_exposures > 1 or self.modifier_k is not None
        if multi and self.family != "gaussian":
            raise ValueError("several exposures or a modifier require family='gaussian'")
        if multi and self.sampler == "gibbs":
            raise ValueError("sampler='gibbs' supports a single exposure without modifier")
        if self.family not in FAMILIES:
            raise ValueError(f"family must be one of {FAMILIES}")
        if self.sampler not in SAMPLERS:
            raise ValueError(f"sampler must be one of {SAMPLERS}")
        if self.sampler == "gibbs" and (self.family != "gaussian" or self.ar1 or self.tau_model != "rw"):
            raise ValueError("sampler='gibbs' supports the Gaussian family without AR(1) and tau_model='rw'")
        if self.ar1 and self.family != "gaussian":
            raise ValueError("ar1=True requires family='gaussian'")
        if self.sampler == "laplace" and self.family != "gaussian":
            raise ValueError("sampler='laplace' requires family='gaussian'")
        if self.parameterization not in ("centred", "noncentred"):
            raise ValueError("parameterization must be 'centred' or 'noncentred'")
        self.structs_ = [build_structure(p, kk, end_constraint=self.end_constraint, lag_scale=self.lag_scale)
                         for p, kk in zip(nlags, ks)]
        self.struct: Structure = self.structs_[0]

    @property
    def nlag_max(self):
        return max(self.nlags_)

    # ----------------------------------------------------------------- fitting
    def fit(self, x, y, covariates=None, offset=None, groups=None, weights=None, modifier=None, *,
            samples=2000, warmup=1000, chains=4, seed=0, thin=1, chain_method="vectorized",
            target_accept=0.9, progress=False, init="laplace", clear_jax_cache=True, dense_mass=True,
            max_tree_depth=10):
        """Fit the model.

        ``x`` is one exposure series, or an array with one column per exposure
        when ``nlag`` is a sequence. ``samples`` are retained per chain after
        ``warmup`` iterations. For ``sampler='nuts'`` the chain is initialised
        at the Laplace mode when ``init='laplace'``. ``groups`` marks separate
        series within ``x`` (for example trading days or assets), so that lags
        never cross a group boundary. ``weights`` are positive observation
        weights (Gaussian family), for example inverse ex-ante variances.
        ``modifier`` is the modifier series (or ``'time'``) when
        ``modifier_k`` is set. Each fit compiles its own JAX executables; with
        ``clear_jax_cache=True`` they are released afterwards so that long
        simulation loops do not exhaust the compiler's memory.
        """
        from .numpyro_model import tau_basis_matrix

        if weights is not None and self.family != "gaussian":
            raise ValueError("weights require family='gaussian'")
        x_arr = np.asarray(x, float)
        if x_arr.ndim == 1:
            x_arr = x_arr[:, None]
        if x_arr.shape[1] != self.n_exposures:
            raise ValueError(f"x must have {self.n_exposures} column(s), one per exposure")
        mod_basis = None
        if self.modifier_k is not None:
            if modifier is None:
                raise ValueError("modifier_k is set; pass modifier= (a series or 'time') to fit")
            kind = "time" if isinstance(modifier, str) else "covariate"
            mod_basis = build_modifier(None if kind == "time" else modifier, K=self.modifier_k, kind=kind)
        self.components_ = [LagComponent(st, mod_basis) for st in self.structs_]
        data, std = build_multi_design(x_arr, y, self.components_, covariates=covariates,
                                       intercept=self.intercept, offset=offset,
                                       standardize=self.standardize, family=self.family,
                                       groups=groups, weights=weights, modifier=modifier)
        self.data_: DesignData = data
        self.std_: Standardizer = std
        self.modifier_ = mod_basis
        # average modifier basis over the training rows: the curve 'at the average state'
        self.mod_weights_ = None if mod_basis is None else mod_basis.basis(data.z).mean(axis=0)
        self.C_ = tau_basis_matrix(self.struct.K_tau, self.priors.tau_basis) if self.tau_model == "basis" else None
        self.fit_info_ = {"sampler": self.sampler, "samples": samples, "warmup": warmup, "chains": chains}
        self._nuts_opts = {"dense_mass": dense_mass, "max_tree_depth": max_tree_depth}

        if self.sampler == "gibbs":
            from .gibbs import run_gibbs

            draws = run_gibbs(data, self.struct, self.priors, samples=samples, warmup=warmup,
                              chains=chains, seed=seed, thin=thin)
            self.draws_ = draws
            self.extra_ = {}
        elif self.family == "gaussian":
            self._fit_gaussian_collapsed(samples, warmup, chains, seed, thin, chain_method,
                                         target_accept, progress, init)
        else:
            self._fit_counts(samples, warmup, chains, seed, thin, chain_method, target_accept, progress)
        self._postprocess()
        if clear_jax_cache and self.sampler != "gibbs":
            import jax

            self.mcmc_ = None
            jax.clear_caches()
        return self

    def _model_and_args(self):
        from .numpyro_model import count_model, gaussian_model

        if self.family == "gaussian":
            args = (self.data_.stats, self.components_, self.data_.q, self.ar1, self.tau_model, self.priors, self.C_)
            return gaussian_model, args
        args = (self.data_.L, self.data_.W, self.data_.y, self.data_.offset, self.struct, self.family,
                self.tau_model, self.priors, self.C_, self.parameterization)
        return count_model, args

    def _fit_gaussian_collapsed(self, samples, warmup, chains, seed, thin, chain_method,
                                target_accept, progress, init):
        from .laplace import laplace_fit
        from .numpyro_model import draw_gaussian_coefficients, ncp_init_values, run_nuts

        model, args = self._model_and_args()
        init_values = None
        if self.sampler == "laplace" or init == "laplace":
            start = ncp_init_values(self.components_, self.tau_model, self.priors, self.C_, self.ar1)
            lap_draws, lap_info = laplace_fit(model, args,
                                              n_draws=samples if self.sampler == "laplace" else 2,
                                              seed=seed, chains=chains if self.sampler == "laplace" else 1,
                                              init_values=start)
            self.laplace_ = lap_info
            mode = lap_info["postprocess"](lap_info["unravel"](lap_info["mode"]))
            self.laplace_["mode_constrained"] = {k: np.asarray(v) for k, v in mode.items()}
            if self.sampler == "nuts":
                init_values = {k: np.asarray(v) for k, v in mode.items() if k in lap_info["latent_names"]}
        if self.sampler == "nuts":
            draws, extra, mcmc = run_nuts(model, args, samples=samples, warmup=warmup, chains=chains,
                                          seed=seed, chain_method=chain_method, target_accept=target_accept,
                                          init_values=init_values, progress=progress, thin=thin,
                                          **self._nuts_opts)
            self.mcmc_ = mcmc
        else:
            draws, extra = lap_draws, {}
            # The Gaussian approximation lives in the non-centred space; guard its
            # tails so that exp(tau) stays representable.
            draws["tau"] = np.clip(draws["tau"], -20.0, 40.0)
        theta, ed = draw_gaussian_coefficients(draws, self.data_.stats, self.components_, self.data_.q,
                                               self.ar1, seed=seed + 1)
        draws["theta"] = theta
        draws["ed"] = ed
        self.draws_, self.extra_ = draws, extra

    def _fit_counts(self, samples, warmup, chains, seed, thin, chain_method, target_accept, progress):
        from .numpyro_model import run_nuts

        model, args = self._model_and_args()
        draws, extra, mcmc = run_nuts(model, args, samples=samples, warmup=warmup, chains=chains,
                                      seed=seed, chain_method=chain_method, target_accept=target_accept,
                                      progress=progress, thin=thin, **self._nuts_opts)
        self.mcmc_ = mcmc
        self.draws_, self.extra_ = draws, extra

    # --------------------------------------------------------------- results
    def _postprocess(self):
        std = self.std_
        theta = self.draws_["theta"]
        x_mean = np.atleast_1d(std.x_mean)
        x_sd = np.atleast_1d(std.x_sd)
        y_scale = std.y_sd if self.family == "gaussian" else 1.0
        self.lag_curves_ = []
        shift = 0.0
        off = 0
        for e, comp in enumerate(self.components_):
            b = theta[..., off: off + comp.dim]
            off += comp.dim
            beta_s = comp.curve(b, weights=self.mod_weights_) if comp.modifier is not None else comp.curve(b)
            sfx = "" if e == 0 else f"_{e}"
            self.draws_["beta_std" + sfx] = beta_s
            self.draws_["lag_curve" + sfx] = beta_s * (y_scale / x_sd[e])
            self.lag_curves_.append(self.draws_["lag_curve" + sfx])
            shift = shift + beta_s.sum(axis=-1) * x_mean[e] / x_sd[e]
        self.n_penalised_ = off
        if self.intercept:
            alpha_s = theta[..., off]
            if self.family == "gaussian":
                self.draws_["intercept"] = std.y_mean + std.y_sd * (alpha_s - shift)
            else:
                self.draws_["intercept"] = alpha_s - shift
        if self.family == "gaussian":
            self.draws_["sigma"] = np.sqrt(self.draws_["sigma2"]) * std.y_sd
        self.lag_curve_ = self.lag_curves_[0]
        self.summary_ = self.lag_curve_summary()
        self.diagnostics_ = self.diagnostics()

    def lag_surface(self, z, exposure=0, probs=(0.025, 0.975)):
        """Posterior summary of a modifier-dependent curve at modifier values
        ``z``: dict with ``z``, ``median``, ``lower``, ``upper``, ``mean`` of shape
        ``(len(z), p + 1)`` on the original scale."""
        comp = self.components_[exposure]
        if comp.modifier is None:
            raise ValueError("this model has no modifier")
        off = sum(c.dim for c in self.components_[:exposure])
        b = self.draws_["theta"][..., off: off + comp.dim].reshape(-1, comp.dim)
        z = np.atleast_1d(np.asarray(z, float))
        surf = comp.curve(b, z=z)                                    # (draws, nz, p + 1)
        scale = (self.std_.y_sd if self.family == "gaussian" else 1.0) / np.atleast_1d(self.std_.x_sd)[exposure]
        surf = surf * scale
        qs = np.quantile(surf, [probs[0], 0.5, probs[1]], axis=0)
        return {"z": z, "lower": qs[0], "median": qs[1], "upper": qs[2], "mean": surf.mean(axis=0),
                "draws": surf}

    @property
    def lags(self):
        return np.arange(self.nlags_[0] + 1)

    def lag_curve_summary(self, probs=(0.025, 0.975), exposure=0):
        p = self.nlags_[exposure]
        lc = self.lag_curves_[exposure].reshape(-1, p + 1)
        qs = np.quantile(lc, [probs[0], 0.5, probs[1]], axis=0)
        return LagCurveSummary(lag=np.arange(p + 1), median=qs[1], lower=qs[0], upper=qs[2],
                               mean=lc.mean(axis=0), probs=tuple(probs))

    def cumulative_effect(self, probs=(0.025, 0.975)):
        """Posterior draws and quantiles of the total effect ``sum_j beta_j``."""
        total = self.lag_curve_.sum(axis=-1).ravel()
        return {"draws": total, "mean": total.mean(), "median": np.median(total),
                "lower": np.quantile(total, probs[0]), "upper": np.quantile(total, probs[1])}

    def peak_lag(self):
        """Posterior distribution of the lag with the largest absolute effect."""
        peaks = np.abs(self.lag_curve_.reshape(-1, self.nlags_[0] + 1)).argmax(axis=1)
        vals, counts = np.unique(peaks, return_counts=True)
        return {"draws": peaks, "mode": int(vals[counts.argmax()]),
                "probabilities": dict(zip(vals.tolist(), (counts / counts.sum()).tolist()))}

    def effective_lag_length(self, eps=None, prob=0.95):
        """Smallest lag beyond which the curve stays within ``eps`` of zero.

        ``eps`` defaults to 5% of the largest absolute posterior-mean coefficient.
        Returns the per-draw lengths and the smallest lag that is an upper
        bound with posterior probability ``prob``.
        """
        p = self.nlags_[0]
        lc = self.lag_curve_.reshape(-1, p + 1)
        if eps is None:
            eps = 0.05 * np.abs(lc.mean(axis=0)).max()
        small = np.abs(lc) < eps
        rev = small[:, ::-1]
        trailing = np.where(rev.all(axis=1), p + 1, np.argmin(rev, axis=1))
        lengths = p + 1 - trailing
        return {"draws": lengths, "eps": float(eps), "length": int(np.quantile(lengths, prob)),
                "median": float(np.median(lengths))}

    def effective_df(self):
        """Posterior draws of the effective degrees of freedom of the lag curve
        (Gaussian family; ``None`` otherwise).

        Unpenalised terms (intercept and covariates) each contribute exactly one
        degree of freedom to the fit and are excluded, so the value describes
        the lag curve alone.
        """
        ed = self.draws_.get("ed")
        if ed is None:
            return None
        ed = np.asarray(ed).ravel() - self.data_.q
        return {"draws": ed, "mean": float(ed.mean()), "sd": float(ed.std())}

    def diagnostics(self):
        """R-hat, effective sample size and divergence counts for the main sites."""
        base = ("tau", "sigma2", "zeta", "phi", "tau_level", "nb_conc", "hs_global", "mod_kappa")
        names = [n for n in self.draws_ if n in base or
                 any(n.startswith(b + "_") and n[len(b) + 1:].isdigit() for b in base)]
        out = summarise_draws(self.draws_, names)
        for e, lc in enumerate(self.lag_curves_):
            key = "lag_curve" if e == 0 else f"lag_curve_{e}"
            out[key] = summarise_draws({key: lc})[key]
        out["cumulative_effect"] = summarise_draws({"c": self.lag_curve_.sum(axis=-1)})["c"]
        if "diverging" in self.extra_:
            out["divergences"] = int(np.asarray(self.extra_["diverging"]).sum())
        return out

    def summary(self, probs=(0.025, 0.975)):
        """Printable summary of the fit."""
        s = self.lag_curve_summary(probs)
        ce = self.cumulative_effect(probs)
        d = self.diagnostics_
        lines = [f"badlm: {self.family} family, tau_model='{self.tau_model}', end_constraint='{self.end_constraint}', "
                 f"sampler='{self.sampler}'",
                 f"nlag={self.nlags_ if self.n_exposures > 1 else self.nlag}, k={self.ks_ if self.n_exposures > 1 else self.k}, "
                 f"chains={self.fit_info_['chains']}, draws per chain={self.fit_info_['samples']}",
                 "",
                 f"cumulative effect: mean {ce['mean']:.4g}, {100*(probs[1]-probs[0]):.0f}% interval "
                 f"[{ce['lower']:.4g}, {ce['upper']:.4g}]"]
        ed = self.effective_df()
        if ed is not None:
            lines.append(f"effective degrees of freedom: {ed['mean']:.2f} (sd {ed['sd']:.2f})")
        pk = self.peak_lag()
        lines.append(f"peak lag (posterior mode): {pk['mode']}")
        el = self.effective_lag_length()
        lines.append(f"effective lag length (95% upper bound, eps={el['eps']:.3g}): {el['length']}")
        if "divergences" in d:
            lines.append(f"divergent transitions: {d['divergences']}")
        worst = max(np.nanmax(v["rhat"]) for k, v in d.items() if isinstance(v, dict))
        least = min(np.nanmin(v["ess"]) for k, v in d.items() if isinstance(v, dict))
        lines.append(f"max R-hat: {worst:.3f}, min bulk ESS: {least:.0f}")
        lines += ["", "lag curve (original scale):", repr(s)]
        return "\n".join(lines)

    # --------------------------------------------------------------- predict
    def predict(self, x_new, covariates=None, offset=None, modifier=None, kind="mean",
                probs=(0.025, 0.975), seed=0):
        """Posterior predictive summaries for a new exposure series (or array
        with one column per exposure).

        Returns a dict with ``mean``, ``lower``, ``upper`` and ``draws`` (shape
        ``(n_draws, len(x_new) - nlag_max)``). ``kind='response'`` adds
        observation noise; ``kind='mean'`` returns the expected response. For a
        modifier model pass ``modifier`` (same length as ``x_new``, or the
        retained length); ``'time'`` uses the end of the training period.
        """
        std = self.std_
        x_new = np.asarray(x_new, float)
        if x_new.ndim == 1:
            x_new = x_new[:, None]
        xs = std.x(x_new)
        p_max = self.nlag_max
        n = x_new.shape[0] - p_max
        z = None
        if self.modifier_ is not None:
            if modifier is None or (isinstance(modifier, str) and modifier == "time"):
                z = np.full(n, self.modifier_.domain[1])
            else:
                z = np.asarray(modifier, float).ravel()
                z = z[p_max:] if z.size == x_new.shape[0] else z
        blocks = []
        for e, comp in enumerate(self.components_):
            p = comp.struct.p
            L = lag_matrix(xs[:, e], p)[p_max - p:]
            blocks.append(comp.design(L, z))
        cols = []
        if self.intercept:
            cols.append(np.ones((n, 1)))
        if covariates is not None:
            Wc = np.asarray(covariates, float).reshape(-1, 1) if np.ndim(covariates) == 1 else np.asarray(covariates, float)
            Wc = Wc[p_max:] if Wc.shape[0] == x_new.shape[0] else Wc
            cols.append(Wc)
        X = np.hstack(blocks + cols) if cols else np.hstack(blocks)
        theta = self.draws_["theta"].reshape(-1, self.draws_["theta"].shape[-1])
        eta = theta @ X.T
        if offset is not None:
            off = np.asarray(offset, float).ravel()
            eta = eta + (off[p_max:] if off.size == x_new.shape[0] else off)
        rng = np.random.default_rng(seed)
        if self.family == "gaussian":
            draws = std.y_mean + std.y_sd * eta
            if kind == "response":
                sigma = np.sqrt(self.draws_["sigma2"]).ravel()[:, None] * std.y_sd
                draws = draws + rng.normal(size=draws.shape) * sigma
        else:
            draws = np.exp(eta)
            if kind == "response":
                if self.family == "poisson":
                    draws = rng.poisson(draws).astype(float)
                else:
                    conc = self.draws_["nb_conc"].ravel()[:, None]
                    draws = rng.negative_binomial(conc, conc / (conc + draws)).astype(float)
        return {"mean": draws.mean(axis=0), "lower": np.quantile(draws, probs[0], axis=0),
                "upper": np.quantile(draws, probs[1], axis=0), "draws": draws}

    # ------------------------------------------------------------------ misc
    def plot(self, ax=None, truth=None, exposure=0, **kwargs):
        from .plotting import plot_lag_curve

        summary = self.summary_ if exposure == 0 else self.lag_curve_summary(exposure=exposure)
        return plot_lag_curve(summary, ax=ax, truth=truth, **kwargs)

    def to_arviz(self):
        """Convert the draws to an ``arviz.InferenceData`` object."""
        import inspect

        import arviz as az

        posterior = {k: np.asarray(v) for k, v in self.draws_.items() if np.asarray(v).ndim >= 2}
        groups = {"posterior": posterior}
        if self.extra_:
            groups["sample_stats"] = {k: np.asarray(v) for k, v in self.extra_.items()}
        if "posterior" in inspect.signature(az.from_dict).parameters:   # ArviZ 0.x: keyword groups
            return az.from_dict(**groups)
        return az.from_dict(groups)                                      # ArviZ >= 1.0: dict of groups
