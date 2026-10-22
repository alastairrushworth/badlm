"""NumPyro model definitions and the NUTS driver."""
from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from numpyro.distributions import constraints
from numpyro.infer import MCMC, NUTS, init_to_value

from .collapsed import (collapsed_loglik, conditional_coefficients, multi_collapsed_loglik,
                        multi_conditional_coefficients)
from .design import bspline_basis
from .priors import Priors

TAU_MODELS = ("rw", "basis", "monotone", "horseshoe", "constant")


def tau_basis_matrix(kt, nb, deg=3):
    """B-spline basis over the knot-gap index used by ``tau_model='basis'``."""
    nb = max(int(nb), deg + 1)
    return bspline_basis(np.arange(kt, dtype=float), nseg=nb - deg, deg=deg, domain=(0.0, kt - 1.0))


def as_components(structs):
    """Normalise ``Structure`` / ``LagComponent`` / list thereof to a list of components."""
    from .design import LagComponent, Structure

    if isinstance(structs, Structure):
        return [LagComponent(structs)]
    if isinstance(structs, LagComponent):
        return [structs]
    return [LagComponent(c) if isinstance(c, Structure) else c for c in structs]


def suffix(e):
    """Site-name suffix of exposure ``e``: the first keeps the original names."""
    return "" if e == 0 else f"_{e}"


def sample_level(pri: Priors, sfx=""):
    """Level of the log precisions: a penalised-complexity prior on the implied
    standard deviation ``exp(-level / 2)`` of coefficient differences."""
    sd = numpyro.sample("tau_level_sd" + sfx, dist.Exponential(pri.tau_level_rate))
    return numpyro.deterministic("tau_level" + sfx, -2.0 * jnp.log(sd))


def sample_tau(kt, tau_model, pri: Priors, C=None, sfx=""):
    """Sample the log precisions ``tau`` (length ``kt``) under one of the prior
    families. Random-walk increments use a non-centred parameterisation."""
    if tau_model == "constant":
        level = sample_level(pri, sfx)
        tau = jnp.full((kt,), level)
    elif tau_model in ("rw", "basis"):
        n_g = kt if tau_model == "rw" else C.shape[1]
        level = sample_level(pri, sfx)
        zeta = numpyro.sample("zeta" + sfx, dist.Exponential(pri.zeta_rate))
        inc = numpyro.sample("tau_inc" + sfx, dist.Normal(0.0, 1.0).expand([n_g - 1]))
        walk = jnp.concatenate([jnp.zeros(1), jnp.cumsum(inc)])
        g = level + zeta * (walk - walk.mean())
        tau = g if tau_model == "rw" else jnp.asarray(C) @ g
    elif tau_model == "monotone":
        level = sample_level(pri, sfx)
        zeta = numpyro.sample("zeta" + sfx, dist.Exponential(pri.zeta_rate))
        e = numpyro.sample("tau_inc" + sfx, dist.Exponential(1.0).expand([kt - 1]))
        tau = level + jnp.concatenate([jnp.zeros(1), jnp.cumsum(zeta * e)])
    elif tau_model == "horseshoe":
        s = numpyro.sample("hs_global" + sfx, dist.HalfCauchy(pri.hs_scale))
        psi = numpyro.sample("hs_local" + sfx, dist.HalfCauchy(1.0).expand([kt]))
        tau = -2.0 * jnp.log(psi * s)
    else:
        raise ValueError(f"tau_model must be one of {TAU_MODELS}")
    return numpyro.deterministic("tau" + sfx, tau)


def sample_kappa(pri: Priors, sfx=""):
    """Precision of coefficient differences along the modifier, through a
    penalised-complexity prior on the corresponding standard deviation."""
    sd = numpyro.sample("mod_sd" + sfx, dist.Exponential(pri.mod_rate))
    return numpyro.deterministic("mod_kappa" + sfx, sd ** -2)


def _init_one(kt, tau_model, pri, C, sfx):
    v = {}
    level_sd = float(np.exp(-0.5 * pri.tau_level_start))
    if tau_model == "constant":
        v["tau_level_sd" + sfx] = level_sd
    elif tau_model in ("rw", "basis"):
        n_g = kt if tau_model == "rw" else C.shape[1]
        v.update({"tau_level_sd" + sfx: level_sd, "zeta" + sfx: 0.5, "tau_inc" + sfx: np.zeros(n_g - 1)})
    elif tau_model == "monotone":
        v.update({"tau_level_sd" + sfx: level_sd, "zeta" + sfx: 0.5, "tau_inc" + sfx: np.full(kt - 1, 0.1)})
    elif tau_model == "horseshoe":
        v.update({"hs_global" + sfx: 0.05, "hs_local" + sfx: np.ones(kt)})
    return v


def ncp_init_values(kt, tau_model, pri: Priors, C=None, ar1=False):
    """A sensible starting point in the model's own (non-centred) sites: a
    constant log precision at the prior mean and a small random-walk scale.
    ``kt`` is the number of log precisions of a single curve, or a list of
    components."""
    v = {"sigma2": 0.5}
    if ar1:
        v["phi"] = 0.0
    if isinstance(kt, (int, np.integer)):
        v.update(_init_one(int(kt), tau_model, pri, C, ""))
        return v
    for e, comp in enumerate(as_components(kt)):
        v.update(_init_one(comp.n_tau, tau_model, pri, C, suffix(e)))
        if comp.modifier is not None:
            v["mod_sd" + suffix(e)] = 0.3
    return v


def gaussian_model(stats, structs, q, ar1, tau_model, pri: Priors, C=None):
    """Collapsed Gaussian model for one or several exposures (``structs`` is a
    ``Structure`` or a list of ``LagComponent``)."""
    comps = as_components(structs)
    taus, kappas = [], []
    for e, comp in enumerate(comps):
        taus.append(sample_tau(comp.n_tau, tau_model, pri, C, suffix(e)))
        kappas.append(sample_kappa(pri, suffix(e)) if comp.modifier is not None else None)
    sigma2 = numpyro.sample("sigma2", dist.InverseGamma(pri.sigma2_a, pri.sigma2_b))
    phi = numpyro.sample("phi", dist.Uniform(-1.0, 1.0)) if ar1 else 0.0
    if len(comps) == 1 and comps[0].modifier is None:
        ll = collapsed_loglik(taus[0], sigma2, phi, stats, comps[0].struct, q, ar1)
    else:
        ll = multi_collapsed_loglik(comps, taus, kappas, sigma2, phi, stats, q, ar1)
    numpyro.factor("loglik", ll / pri.temper)


def count_model(L, W, y, offset, struct, family, tau_model, pri: Priors, C=None,
                parameterization="centred"):
    """Poisson or negative binomial response.

    ``parameterization='centred'`` samples the reduced coefficients directly
    under their Gaussian prior, which mixes well when the data are informative
    about the curve; ``'noncentred'`` samples standard normal increments
    ``w`` with ``d_k = exp(-tau_k / 2) w_k``, which is preferable when they are
    not.
    """
    kt = struct.K_tau
    tau = sample_tau(kt, tau_model, pri, C)
    D_r = jnp.asarray(struct.D_r)
    Z = jnp.asarray(struct.Z)
    if parameterization == "centred":
        if struct.m >= 1:
            Q_r = (D_r.T * jnp.exp(tau)) @ D_r
            b_r = numpyro.sample("b_r", dist.MultivariateNormal(jnp.zeros(struct.K_r), precision_matrix=Q_r))
        else:
            b_r = numpyro.sample("b_r", dist.ImproperUniform(constraints.real, (), (struct.K_r,)))
            numpyro.factor("b_prior", 0.5 * (struct.logdet_const + jnp.sum(tau))
                           - 0.5 * jnp.sum(jnp.exp(tau) * (D_r @ b_r) ** 2))
        b = Z @ b_r
    elif parameterization == "noncentred":
        w = numpyro.sample("w", dist.Normal(0.0, 1.0).expand([kt]))
        scale = jnp.exp(-0.5 * tau)
        if struct.m == 2:
            a = jnp.asarray(struct.a_slope) * scale
            w = w - a * jnp.dot(a, w) / jnp.dot(a, a)
        d = scale * w
        b_raw = jnp.concatenate([jnp.zeros(1), jnp.cumsum(d)])
        if struct.m == 0:
            level = numpyro.sample("level", dist.Normal(0.0, pri.coef_sd))
        else:
            level = -jnp.dot(jnp.asarray(struct.C[0]), b_raw)
        b = level + b_raw
        b_r = Z.T @ b
    else:
        raise ValueError("parameterization must be 'centred' or 'noncentred'")
    beta = numpyro.deterministic("beta", jnp.asarray(struct.B) @ b)
    eta = jnp.asarray(L) @ beta + jnp.asarray(offset)
    q = W.shape[1]
    if q:
        gamma = numpyro.sample("gamma", dist.Normal(0.0, pri.coef_sd).expand([q]))
        eta = eta + jnp.asarray(W) @ gamma
        numpyro.deterministic("theta", jnp.concatenate([b_r, gamma]))
    else:
        numpyro.deterministic("theta", b_r)
    mean = jnp.exp(eta)
    if family == "poisson":
        numpyro.sample("y", dist.Poisson(mean), obs=jnp.asarray(y))
    elif family == "negbin":
        conc = numpyro.sample("nb_conc", dist.HalfCauchy(5.0))
        numpyro.sample("y", dist.NegativeBinomial2(mean, conc), obs=jnp.asarray(y))
    else:
        raise ValueError("family must be 'gaussian', 'poisson' or 'negbin'")


def run_nuts(model, model_args, *, samples=2000, warmup=1000, chains=4, seed=0,
             chain_method="vectorized", target_accept=0.9, init_values=None, progress=False,
             thin=1, dense_mass=True, max_tree_depth=10):
    opts = dict(target_accept_prob=target_accept, dense_mass=dense_mass, max_tree_depth=max_tree_depth)
    if init_values:
        kernel = NUTS(model, init_strategy=init_to_value(values=init_values), **opts)
    else:
        kernel = NUTS(model, **opts)
    mcmc = MCMC(kernel, num_warmup=warmup, num_samples=samples, num_chains=chains,
                chain_method=chain_method, progress_bar=progress, thinning=thin)
    mcmc.run(jax.random.PRNGKey(seed), *model_args, extra_fields=("diverging", "num_steps"))
    draws = {k: np.asarray(v) for k, v in mcmc.get_samples(group_by_chain=True).items()}
    extra = {k: np.asarray(v) for k, v in mcmc.get_extra_fields(group_by_chain=True).items()}
    return draws, extra, mcmc


def hyper_draws(draws, structs):
    """Per-component ``tau`` and ``kappa`` arrays from a draws dictionary."""
    comps = as_components(structs)
    taus = [jnp.asarray(draws["tau" + suffix(e)]) for e in range(len(comps))]
    kappas = [jnp.asarray(draws["mod_kappa" + suffix(e)]) if c.modifier is not None else None
              for e, c in enumerate(comps)]
    return comps, taus, kappas


def draw_gaussian_coefficients(draws, stats, structs, q, ar1, seed=0):
    """Post-hoc draws of the coefficients and effective degrees of freedom for
    every retained hyperparameter draw of a collapsed Gaussian fit."""
    comps, taus, kappas = hyper_draws(draws, structs)
    sigma2 = jnp.asarray(draws["sigma2"])
    phi = jnp.asarray(draws["phi"]) if ar1 else jnp.zeros_like(sigma2)
    chains, n = sigma2.shape
    dim = sum(c.dim for c in comps) + q
    z = jax.random.normal(jax.random.PRNGKey(seed), (chains, n, dim))
    has_mod = [c.modifier is not None for c in comps]
    kap_in = [k if h else jnp.zeros_like(sigma2) for k, h in zip(kappas, has_mod)]

    def one(ts, ks, s2, ph, zz):
        ks = [k if h else None for k, h in zip(ks, has_mod)]
        if len(comps) == 1 and not has_mod[0]:
            return conditional_coefficients(ts[0], s2, ph, stats, comps[0].struct, q, zz)
        return multi_conditional_coefficients(comps, ts, ks, s2, ph, stats, q, zz)

    f = jax.jit(jax.vmap(jax.vmap(one)))
    theta, ed = f(taus, kap_in, sigma2, phi, z)
    return np.asarray(theta), np.asarray(ed)
