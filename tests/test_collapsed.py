import jax.numpy as jnp
import numpy as np
import pytest
from scipy.stats import multivariate_normal

from badlm.collapsed import collapsed_loglik, conditional_coefficients
from badlm.design import ar1_whitening_matrix, build_design, build_structure


@pytest.fixture
def small(rng):
    st = build_structure(15, 8, end_constraint="zero")
    x = rng.normal(size=120)
    y = rng.normal(size=120)
    data, _ = build_design(x, y, st, intercept=False, standardize=False)
    return st, data


def test_collapsed_matches_brute_force_marginal(small, rng):
    st, data = small
    for _ in range(3):
        tau = rng.normal(3, 1, st.K_tau)
        s2 = rng.uniform(0.5, 2)
        cov = s2 * np.eye(data.N) + data.X @ np.linalg.solve(st.lam_matrix(tau), data.X.T)
        brute = multivariate_normal(np.zeros(data.N), cov).logpdf(data.y)
        ours = float(collapsed_loglik(jnp.asarray(tau), s2, 0.0, data.stats, st, 0, False))
        assert abs(brute - ours) < 1e-6


def test_collapsed_matches_brute_force_with_ar1(small, rng):
    st, data = small
    phi, s2 = 0.37, 1.3
    L = ar1_whitening_matrix(data.N, phi)
    tau = rng.normal(3, 1, st.K_tau)
    cov = s2 * np.linalg.inv(L.T @ L) + data.X @ np.linalg.solve(st.lam_matrix(tau), data.X.T)
    brute = multivariate_normal(np.zeros(data.N), cov).logpdf(data.y)
    ours = float(collapsed_loglik(jnp.asarray(tau), s2, phi, data.stats, st, 0, True))
    assert abs(brute - ours) < 1e-6


def test_flat_prior_intercept_is_limit_of_vague_prior(rng):
    st = build_structure(15, 8, end_constraint="zero")
    x, y = rng.normal(size=120), rng.normal(size=120)
    data, _ = build_design(x, y, st, intercept=True, standardize=False)
    tau_a = rng.normal(3, 1, st.K_tau)
    tau_b = tau_a + rng.normal(0, 0.5, st.K_tau)
    s2 = 1.3

    def brute(tau, v=1e7):
        Qf = np.zeros((st.K_r + 1, st.K_r + 1))
        Qf[:st.K_r, :st.K_r] = st.lam_matrix(tau)
        Qf[-1, -1] = 1 / v
        cov = s2 * np.eye(data.N) + data.X @ np.linalg.solve(Qf, data.X.T)
        return multivariate_normal(np.zeros(data.N), cov).logpdf(data.y)

    ours = [float(collapsed_loglik(jnp.asarray(t), s2, 0.0, data.stats, st, 1, False)) for t in (tau_a, tau_b)]
    assert abs((brute(tau_a) - brute(tau_b)) - (ours[0] - ours[1])) < 1e-4


def test_conditional_coefficients_mean_and_ed(small, rng):
    st, data = small
    tau = rng.normal(3, 1, st.K_tau)
    s2 = 0.8
    Q = st.lam_matrix(tau)
    A = Q + data.stats.S_a / s2
    mu = np.linalg.solve(A, data.stats.s_a / s2)
    theta, ed = conditional_coefficients(jnp.asarray(tau), s2, 0.0, data.stats, st, 0, jnp.zeros(st.K_r))
    assert np.allclose(np.asarray(theta), mu)
    ed_ref = np.trace(np.linalg.solve(data.stats.S_a + s2 * Q, data.stats.S_a))
    assert np.isclose(float(ed), ed_ref)
    assert 0 < float(ed) <= st.K_r
