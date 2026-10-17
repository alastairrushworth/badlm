# badlm

Bayesian adaptive distributed lag models in Python.

`badlm` infers the temporally delayed dependence between an exposure series
and a response. It implements the model of

- Rushworth, A. (2018). *Bayesian adaptive distributed lag models.*
  [arXiv:1801.06670](https://arxiv.org/abs/1801.06670). A revised copy of the
  paper source lives in [`paper/`](paper/).

The lag curve is represented by a B-spline basis whose coefficients receive a
first-difference penalty with **lag-varying precision**; the log precisions are
themselves smoothed, so the amount of smoothing adapts to the shape of the curve
and the maximum lag can be chosen generously without overfitting.

This version is a complete rewrite of the original R package. The spline
coefficients are integrated out analytically and the remaining low-dimensional
posterior is sampled with the No-U-Turn Sampler (via NumPyro and JAX), which
removes the funnel between the curve and its smoothing parameters that slowed
the original Gibbs/Metropolis scheme.

## Installation

```bash
pip install git+https://github.com/alastairrushworth/badlm
```

Python 3.10 or later; depends on NumPy, SciPy, JAX, NumPyro and matplotlib.
Optional extras: `pip install "badlm[arviz,pandas]"`.

## Quick start

```python
import numpy as np
from badlm import BADLM, readme_lag_curve, simulate_dlm

# a 'hump' lag curve over 50 lags and data simulated from it
truth = readme_lag_curve(50)
x, y = simulate_dlm(truth, n=500, standardize_x=True, ar_eps=0.0, sd_eps=0.01, seed=2018)

model = BADLM(nlag=50, k=30)               # adaptive model, exact zero constraint at lag 50
model.fit(x, y, samples=1000, warmup=1000, chains=4, seed=1)
print(model.summary())
model.plot(truth=truth)
```

![fitted lag curve](figures/readme_fit.png)

The fitted object exposes the posterior draws of the lag curve on the original
scale (`model.lag_curve_`, shape `(chains, draws, nlag + 1)`) and a set of derived
quantities:

```python
model.lag_curve_summary()        # median and 95% band per lag (to_dataframe() available)
model.cumulative_effect()        # total effect of a unit change in x
model.peak_lag()                 # posterior distribution of the lag with largest effect
model.effective_lag_length()     # lag beyond which the curve is effectively zero
model.effective_df()             # effective degrees of freedom of the fitted curve
model.diagnostics()              # R-hat, bulk ESS, divergent transitions
model.predict(x_new)             # posterior predictive mean / response
model.to_arviz()                 # arviz.InferenceData for further analysis
```

Adaptive versus non-adaptive smoothing on the same data (`tau_model='constant'`
is model M1 of the paper; the default `'rw'` is M3):

![adaptive versus non-adaptive](figures/readme_adaptive_vs_not.png)

## Model options

```python
BADLM(
    nlag=50,                 # maximum lag p
    k=None,                  # B-spline functions; default ceil(2(p + 1)/3)
    tau_model="rw",          # 'rw' | 'basis' | 'monotone' | 'horseshoe' | 'constant'
    end_constraint="zero",   # 'zero' | 'zero_slope' | 'none'
    family="gaussian",       # 'gaussian' | 'poisson' | 'negbin'
    ar1=False,               # AR(1) errors (Gaussian family)
    intercept=True,
    standardize=True,        # internal standardisation; results on the original scale
    sampler="nuts",          # 'nuts' | 'gibbs' | 'laplace'
    priors=Priors(),         # hyperparameters, see badlm.priors.Priors
)
```

**Priors on the log precisions.** `'rw'` is the random walk of the paper with a
weakly informative prior on its level and a penalised-complexity (exponential)
prior on its scale, so the model shrinks towards the non-adaptive fit unless the
data ask for adaptivity. `'basis'` represents the log precisions with a handful
of B-splines, `'monotone'` forces smoothing to increase with lag, `'horseshoe'`
uses independent local shrinkage of each difference, and `'constant'` is the
non-adaptive model.

**End constraints.** The curve can be forced to zero at the maximum lag, or to
zero with zero slope, by projecting the coefficients onto the constraint null
space; no tuning parameter is involved and the closed-form generalised
determinant used in the likelihood is unchanged.

**Samplers.** `'nuts'` integrates the coefficients out and samples the log
precisions, error variance and AR(1) coefficient with NUTS; each likelihood
evaluation costs O(k^3) irrespective of the series length because only
sufficient statistics are used. Chains start at the Laplace mode. `'gibbs'` is a
Gibbs sampler with slice updates and interweaving, kept as an independent check
of the NUTS implementation. `'laplace'` draws from a Gaussian approximation at
the mode and takes a few seconds.

**Other families.** Poisson and negative binomial responses use a non-centred
parameterisation of the coefficients and NUTS on the full posterior; an
`offset` can be supplied to `fit`. Additional linear covariates are passed as
`covariates`.

## Validation

The test suite (`pytest`) checks the collapsed likelihood against brute-force
Gaussian marginals, including AR(1) errors; the closed-form generalised
determinants; a Geweke "getting it right" test of the Gibbs sampler; agreement
between the NUTS, Gibbs and Laplace fits; and recovery of known lag curves for
every prior, constraint and family. `pytest -m slow` additionally runs a
simulation-based calibration check of the NUTS sampler.

## Replicating the paper's comparison

`examples/replicate_paper.py` fits the adaptive (M3) and non-adaptive (M1)
models to simulated data under five lag-curve shapes, as in Section 3.1 of the
paper. The published curves are not available in closed form, so
`badlm.paper_lag_curves` provides approximations with the same features.

Run `python examples/replicate_paper.py 10` to produce a Markdown table of RMSE,
squared bias and effective degrees of freedom for both models.

## Citation

```
@article{rushworth2018badlm,
  title   = {Bayesian adaptive distributed lag models},
  author  = {Rushworth, Alastair},
  journal = {arXiv preprint arXiv:1801.06670},
  year    = {2018}
}
```
