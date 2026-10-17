"""User-facing estimator."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .design import DesignData, Standardizer, Structure, build_design, build_structure, lag_matrix
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
    nlag : int
        Maximum lag ``p``. With adaptive smoothing it can be chosen generously.
    k : int or None
        Number of cubic B-spline functions for the lag curve. Defaults to
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
    """

    nlag: int
    k: int | None = None
    tau_model: str = "rw"
    end_constraint: str = "zero"
    family: str = "gaussian"
    ar1: bool = False
    intercept: bool = True
    standardize: bool = True
    sampler: str = "nuts"
    priors: Priors = field(default_factory=Priors)

    def __post_init__(self):
        if self.k is None:
            self.k = int(np.ceil(2.0 * (self.nlag + 1) / 3.0))
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
        self.struct: Structure = build_structure(self.nlag, self.k, end_constraint=self.end_constraint)

    # ----------------------------------------------------------------- fitting
    def fit(self, x, y, covariates=None, offset=None, *, samples=2000, warmup=1000, chains=4,
            seed=0, thin=1, chain_method="vectorized", target_accept=0.9, progress=False,
            init="laplace", clear_jax_cache=True):
        """Fit the model.

        ``samples`` are retained per chain after ``warmup`` iterations. For
        ``sampler='nuts'`` the chain is initialised at the Laplace mode when
        ``init='laplace'``. Each fit compiles its own JAX executables; with
        ``clear_jax_cache=True`` they are released afterwards so that long
        simulation loops do not exhaust the compiler's memory.
        """
        from .numpyro_model import tau_basis_matrix

        data, std = build_design(x, y, self.struct, covariates=covariates, intercept=self.intercept,
                                 offset=offset, standardize=self.standardize, family=self.family)
        self.data_: DesignData = data
        self.std_: Standardizer = std
        self.C_ = tau_basis_matrix(self.struct.K_tau, self.priors.tau_basis) if self.tau_model == "basis" else None
        self.fit_info_ = {"sampler": self.sampler, "samples": samples, "warmup": warmup, "chains": chains}

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
            args = (self.data_.stats, self.struct, self.data_.q, self.ar1, self.tau_model, self.priors, self.C_)
            return gaussian_model, args
        args = (self.data_.L, self.data_.W, self.data_.y, self.data_.offset, self.struct, self.family,
                self.tau_model, self.priors, self.C_)
        return count_model, args

    def _fit_gaussian_collapsed(self, samples, warmup, chains, seed, thin, chain_method,
                                target_accept, progress, init):
        from .laplace import laplace_fit
        from .numpyro_model import draw_gaussian_coefficients, ncp_init_values, run_nuts

        model, args = self._model_and_args()
        init_values = None
        if self.sampler == "laplace" or init == "laplace":
            start = ncp_init_values(self.struct.K_tau, self.tau_model, self.priors, self.C_, self.ar1)
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
                                          init_values=init_values, progress=progress, thin=thin)
            self.mcmc_ = mcmc
        else:
            draws, extra = lap_draws, {}
            # The Gaussian approximation lives in the non-centred space; guard its
            # tails so that exp(tau) stays representable.
            draws["tau"] = np.clip(draws["tau"], -20.0, 40.0)
        theta, ed = draw_gaussian_coefficients(draws, self.data_.stats, self.struct, self.data_.q,
                                               self.ar1, seed=seed + 1)
        draws["theta"] = theta
        draws["ed"] = ed
        self.draws_, self.extra_ = draws, extra

    def _fit_counts(self, samples, warmup, chains, seed, thin, chain_method, target_accept, progress):
        from .numpyro_model import run_nuts

        model, args = self._model_and_args()
        draws, extra, mcmc = run_nuts(model, args, samples=samples, warmup=warmup, chains=chains,
                                      seed=seed, chain_method=chain_method, target_accept=target_accept,
                                      progress=progress, thin=thin)
        self.mcmc_ = mcmc
        self.draws_, self.extra_ = draws, extra

    # --------------------------------------------------------------- results
    def _postprocess(self):
        st, std = self.struct, self.std_
        theta = self.draws_["theta"]
        beta_s = theta[..., : st.K_r] @ st.B_r.T                       # (chains, draws, p + 1)
        if self.family == "gaussian":
            scale = std.y_sd / std.x_sd
        else:
            scale = 1.0 / std.x_sd
        self.draws_["beta_std"] = beta_s
        self.draws_["lag_curve"] = beta_s * scale
        if self.intercept:
            alpha_s = theta[..., st.K_r]
            shift = beta_s.sum(axis=-1) * std.x_mean / std.x_sd
            if self.family == "gaussian":
                self.draws_["intercept"] = std.y_mean + std.y_sd * (alpha_s - shift)
            else:
                self.draws_["intercept"] = alpha_s - shift
        if self.family == "gaussian":
            self.draws_["sigma"] = np.sqrt(self.draws_["sigma2"]) * std.y_sd
        self.lag_curve_ = self.draws_["lag_curve"]
        self.summary_ = self.lag_curve_summary()
        self.diagnostics_ = self.diagnostics()

    @property
    def lags(self):
        return np.arange(self.nlag + 1)

    def lag_curve_summary(self, probs=(0.025, 0.975)):
        lc = self.lag_curve_.reshape(-1, self.nlag + 1)
        qs = np.quantile(lc, [probs[0], 0.5, probs[1]], axis=0)
        return LagCurveSummary(lag=self.lags, median=qs[1], lower=qs[0], upper=qs[2],
                               mean=lc.mean(axis=0), probs=tuple(probs))

    def cumulative_effect(self, probs=(0.025, 0.975)):
        """Posterior draws and quantiles of the total effect ``sum_j beta_j``."""
        total = self.lag_curve_.sum(axis=-1).ravel()
        return {"draws": total, "mean": total.mean(), "median": np.median(total),
                "lower": np.quantile(total, probs[0]), "upper": np.quantile(total, probs[1])}

    def peak_lag(self):
        """Posterior distribution of the lag with the largest absolute effect."""
        peaks = np.abs(self.lag_curve_.reshape(-1, self.nlag + 1)).argmax(axis=1)
        vals, counts = np.unique(peaks, return_counts=True)
        return {"draws": peaks, "mode": int(vals[counts.argmax()]),
                "probabilities": dict(zip(vals.tolist(), (counts / counts.sum()).tolist()))}

    def effective_lag_length(self, eps=None, prob=0.95):
        """Smallest lag beyond which the curve stays within ``eps`` of zero.

        ``eps`` defaults to 5% of the largest absolute posterior-mean coefficient.
        Returns the per-draw lengths and the smallest lag that is an upper
        bound with posterior probability ``prob``.
        """
        lc = self.lag_curve_.reshape(-1, self.nlag + 1)
        if eps is None:
            eps = 0.05 * np.abs(lc.mean(axis=0)).max()
        small = np.abs(lc) < eps
        rev = small[:, ::-1]
        trailing = np.where(rev.all(axis=1), self.nlag + 1, np.argmin(rev, axis=1))
        lengths = self.nlag + 1 - trailing
        return {"draws": lengths, "eps": float(eps), "length": int(np.quantile(lengths, prob)),
                "median": float(np.median(lengths))}

    def effective_df(self):
        """Posterior draws of the effective degrees of freedom of the lag curve
        (Gaussian family; ``None`` otherwise)."""
        ed = self.draws_.get("ed")
        if ed is None:
            return None
        ed = np.asarray(ed).ravel()
        return {"draws": ed, "mean": float(ed.mean()), "sd": float(ed.std())}

    def diagnostics(self):
        """R-hat, effective sample size and divergence counts for the main sites."""
        names = [n for n in ("tau", "sigma2", "zeta", "phi", "tau_level", "nb_conc", "hs_global") if n in self.draws_]
        out = summarise_draws(self.draws_, names)
        out["lag_curve"] = summarise_draws({"lag_curve": self.lag_curve_})["lag_curve"]
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
                 f"nlag={self.nlag}, k={self.k}, chains={self.fit_info_['chains']}, draws per chain={self.fit_info_['samples']}",
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
    def predict(self, x_new, covariates=None, offset=None, kind="mean", probs=(0.025, 0.975), seed=0):
        """Posterior predictive summaries for a new exposure series.

        Returns a dict with ``mean``, ``lower``, ``upper`` and ``draws`` (shape
        ``(n_draws, len(x_new) - nlag)``). ``kind='response'`` adds observation
        noise; ``kind='mean'`` returns the expected response.
        """
        st, std = self.struct, self.std_
        x_new = np.asarray(x_new, float).ravel()
        L = lag_matrix(std.x(x_new), self.nlag)
        n = L.shape[0]
        theta = self.draws_["theta"].reshape(-1, st.K_r + self.data_.q)
        eta = theta[:, : st.K_r] @ (L @ st.B_r).T
        cols = []
        if self.intercept:
            cols.append(np.ones((n, 1)))
        if covariates is not None:
            Wc = np.asarray(covariates, float).reshape(-1, 1) if np.ndim(covariates) == 1 else np.asarray(covariates, float)
            Wc = Wc[self.nlag:] if Wc.shape[0] == x_new.size else Wc
            cols.append(Wc)
        if cols:
            W = np.hstack(cols)
            eta = eta + theta[:, st.K_r:] @ W.T
        if offset is not None:
            off = np.asarray(offset, float).ravel()
            eta = eta + (off[self.nlag:] if off.size == x_new.size else off)
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
    def plot(self, ax=None, truth=None, **kwargs):
        from .plotting import plot_lag_curve

        return plot_lag_curve(self.summary_, ax=ax, truth=truth, **kwargs)

    def to_arviz(self):
        """Convert the draws to an ``arviz.InferenceData`` object."""
        import arviz as az

        posterior = {k: np.asarray(v) for k, v in self.draws_.items() if np.asarray(v).ndim >= 2}
        groups = {"posterior": posterior}
        if self.extra_:
            groups["sample_stats"] = {k: np.asarray(v) for k, v in self.extra_.items()}
        try:                                    # ArviZ >= 1.0 takes a dict of groups
            return az.from_dict(groups)
        except TypeError:                       # ArviZ 0.x takes keyword groups
            return az.from_dict(**groups)
