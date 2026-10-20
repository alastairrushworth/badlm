"""Maximum a posteriori fit and Laplace approximation of the collapsed posterior."""
from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
from jax.flatten_util import ravel_pytree
from numpyro.infer import init_to_value
from numpyro.infer.util import initialize_model
from scipy.optimize import minimize


def laplace_fit(model, model_args, *, n_draws=2000, seed=0, chains=1, init_values=None):
    """Find the posterior mode in NumPyro's unconstrained space, form a Gaussian
    approximation from the Hessian there, and return draws from it in the same
    ``(chains, draws, ...)`` layout as ``run_nuts``."""
    rng = jax.random.PRNGKey(seed)
    strategy = init_to_value(values=init_values) if init_values else None
    info = (initialize_model(rng, model, model_args=model_args, init_strategy=strategy) if strategy
            else initialize_model(rng, model, model_args=model_args))
    z0, unravel = ravel_pytree(info.param_info.z)
    pot = jax.jit(lambda z: info.potential_fn(unravel(z)))
    grad = jax.jit(jax.grad(pot))
    hess = jax.jit(jax.hessian(pot))

    big = 1e30

    def f(z):
        val = float(pot(jnp.asarray(z)))
        return val if np.isfinite(val) else big

    def g(z):
        gr = np.asarray(grad(jnp.asarray(z)), dtype=float)
        return gr if np.all(np.isfinite(gr)) else np.zeros_like(gr)

    res = minimize(f, np.asarray(z0, float), jac=g, method="L-BFGS-B",
                   options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8})
    mode = res.x
    H = np.asarray(hess(jnp.asarray(mode)), float)
    H = 0.5 * (H + H.T)
    evals, evecs = np.linalg.eigh(H)
    evals = np.maximum(evals, 1e-8 * max(evals.max(), 1e-8))
    cov = (evecs / evals) @ evecs.T
    rng_np = np.random.default_rng(seed)
    z = rng_np.multivariate_normal(mode, cov, size=(chains, n_draws))
    post = jax.jit(jax.vmap(jax.vmap(lambda zz: info.postprocess_fn(unravel(zz)))))
    draws = {k: np.asarray(v) for k, v in post(jnp.asarray(z)).items()}
    draws.pop("loglik", None)
    return draws, {"mode": mode, "cov": cov, "converged": bool(res.success),
                   "neg_log_post": float(res.fun), "unravel": unravel,
                   "postprocess": info.postprocess_fn, "latent_names": list(info.param_info.z.keys())}
