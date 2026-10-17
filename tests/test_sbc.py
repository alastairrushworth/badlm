"""Simulation-based calibration (Talts et al., 2018) of the collapsed NUTS
sampler: for data simulated from the prior, the rank of the true parameter among
posterior draws is uniform. Slow; run with ``pytest -m slow``."""
import numpy as np
import pytest
from scipy.stats import kstest

from badlm import BADLM
from badlm.design import build_structure, lag_matrix
from badlm.priors import Priors, simulate_tau


def _sym(a):
    return 0.5 * (a + a.T)


@pytest.mark.slow
def test_sbc_collapsed_nuts():
    rng = np.random.default_rng(2018)
    p, K = 10, 7
    pri = Priors(sigma2_a=3.0, sigma2_b=2.0, tau_level_mean=1.0, tau_level_sd=1.0,
                 zeta_scale=1.0, zeta_tail_prob=0.3)
    st = build_structure(p, K, end_constraint="zero")
    x = rng.normal(size=120)
    L = lag_matrix(x, p)
    n_rep, n_draw = 40, 200
    ranks = {"sigma2": [], "tau_mean": [], "beta_peak": []}
    for r in range(n_rep):
        tau, zeta = simulate_tau(rng, st.K_tau, pri)
        sigma2 = 1.0 / rng.gamma(pri.sigma2_a, 1.0 / pri.sigma2_b)
        b = rng.multivariate_normal(np.zeros(st.K_r), _sym(np.linalg.inv(st.lam_matrix(tau))))
        beta = st.B_r @ b
        y = L @ beta + np.sqrt(sigma2) * rng.normal(size=L.shape[0])
        m = BADLM(nlag=p, k=K, intercept=False, standardize=False, priors=pri)
        m.fit(x, y, samples=n_draw * 4, warmup=300, chains=1, seed=r, thin=4)
        ranks["sigma2"].append((m.draws_["sigma2"].ravel() < sigma2).sum())
        ranks["tau_mean"].append((m.draws_["tau"].reshape(-1, st.K_tau).mean(1) < tau.mean()).sum())
        ranks["beta_peak"].append((m.lag_curve_.reshape(-1, p + 1)[:, 3] < beta[3]).sum())
    for name, rk in ranks.items():
        u = (np.asarray(rk) + rng.uniform(size=n_rep)) / (n_draw + 1)
        assert kstest(u, "uniform").pvalue > 0.01, (name, rk)
