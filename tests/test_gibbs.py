"""'Getting it right' (Geweke, 2004) validation of the Gibbs sampler.

The marginal-conditional simulator draws (theta, y) from the prior; the
successive-conditional simulator alternates one Gibbs sweep of theta | y with a
draw of y | theta. If the sampler leaves the posterior invariant, both produce
the same joint distribution, so moments of test functions must agree.
"""
import numpy as np
import pytest

from badlm.design import build_design, build_structure
from badlm.diagnostics import ess
from badlm.gibbs import GibbsSampler
from badlm.priors import Priors, simulate_tau


def _sym(a):
    return 0.5 * (a + a.T)


def _prior_draw(rng, st, pri):
    tau, zeta = simulate_tau(rng, st.K_tau, pri)
    sigma2 = 1.0 / rng.gamma(pri.sigma2_a, 1.0 / pri.sigma2_b)
    b = rng.multivariate_normal(np.zeros(st.K_r), _sym(np.linalg.inv(st.lam_matrix(tau))))
    return tau, zeta, sigma2, b


def _funcs(tau, zeta, sigma2, b):
    return np.concatenate([tau, [np.log(zeta), np.log(sigma2), tau.mean() ** 2], b, b**2])


@pytest.mark.parametrize("interweave", [False, True])
def test_getting_it_right(interweave):
    rng = np.random.default_rng(7 + int(interweave))
    st = build_structure(8, 6, end_constraint="zero")          # proper prior on b_r
    pri = Priors(sigma2_a=3.0, sigma2_b=2.0, tau_level_scale=0.5, tau_level_tail_prob=0.2,
                 zeta_scale=1.0, zeta_tail_prob=0.3)
    x = rng.normal(size=30)
    data, _ = build_design(x, rng.normal(size=30), st, intercept=False, standardize=False)
    X = data.X
    M = 3000
    A = np.array([_funcs(*_prior_draw(rng, st, pri)) for _ in range(M)])

    sampler = GibbsSampler(data, st, pri, rng, interweave=interweave)
    tau, zeta, sigma2, b = _prior_draw(rng, st, pri)
    sampler.tau, sampler.zeta, sampler.sigma2, sampler.theta = tau, zeta, sigma2, b
    sampler._cache_prior()
    B = np.empty_like(A)
    for i in range(M):
        y = X @ sampler.theta + np.sqrt(sampler.sigma2) * rng.normal(size=data.N)
        sampler.set_y(y)
        sampler.step()
        B[i] = _funcs(sampler.tau, sampler.zeta, sampler.sigma2, sampler.theta)

    se_a = A.std(axis=0) / np.sqrt(M)
    se_b = B.std(axis=0) / np.sqrt(np.maximum(ess(B[None]), 10.0))
    z = (A.mean(axis=0) - B.mean(axis=0)) / np.sqrt(se_a**2 + se_b**2)
    assert np.abs(z).max() < 4.5, z


def test_gibbs_preserves_end_constraint(rng):
    st = build_structure(10, 7, end_constraint="zero")
    x = rng.normal(size=80)
    data, _ = build_design(x, rng.normal(size=80), st)
    sampler = GibbsSampler(data, st, Priors(), rng)
    for _ in range(20):
        sampler.step()
        b = st.Z @ sampler.theta[:st.K_r]
        assert abs(st.B[-1] @ b) < 1e-10
        assert np.all(np.isfinite(sampler.tau))
