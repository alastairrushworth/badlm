import numpy as np
import pytest

from badlm import BADLM, lag_matrix, readme_lag_curve, simulate_dlm


@pytest.fixture(scope="module")
def hump_data():
    truth = readme_lag_curve(20)
    x, y = simulate_dlm(truth, n=400, standardize_x=True, ar_eps=0.0, sd_eps=0.01, seed=11)
    return truth, x, y


def _check_recovery(m, truth, rmse_frac=0.2, min_cover=0.7):
    s = m.summary_
    rmse = np.sqrt(np.mean((s.median - truth) ** 2))
    assert rmse < rmse_frac * np.abs(truth).max()
    assert np.mean((truth >= s.lower) & (truth <= s.upper)) >= min_cover
    assert np.all(s.lower <= s.median) and np.all(s.median <= s.upper)


def test_nuts_recovers_lag_curve_and_summaries(hump_data):
    truth, x, y = hump_data
    m = BADLM(nlag=20, k=12).fit(x, y, samples=400, warmup=400, chains=2, seed=1)
    _check_recovery(m, truth)
    assert m.lag_curve_.shape == (2, 400, 21)
    d = m.diagnostics_
    assert np.nanmax(d["lag_curve"]["rhat"]) < 1.1
    assert d["divergences"] <= 8                      # 1% of 800 draws
    ce = m.cumulative_effect(probs=(0.005, 0.995))
    assert ce["lower"] <= truth.sum() <= ce["upper"]
    assert abs(ce["mean"] - truth.sum()) < 0.05 * truth.sum()
    assert abs(m.peak_lag()["mode"] - 10) <= 2
    assert 5 <= m.effective_lag_length()["length"] <= 21
    ed = m.effective_df()
    assert 1 < ed["mean"] < 11
    text = m.summary()
    assert "cumulative effect" in text and "lag curve" in text
    df = m.lag_curve_summary().to_dataframe()
    assert list(df.columns)[:2] == ["lag", "mean"]
    pr = m.predict(x[:100])
    assert pr["mean"].shape == (80,) and pr["draws"].shape == (800, 80)
    pr2 = m.predict(x[:100], kind="response")
    assert np.all(pr2["upper"] >= pr2["lower"])
    idata = m.to_arviz()
    assert "lag_curve" in idata.posterior
    ax = m.plot(truth=truth)
    assert ax.get_xlabel() == "lag"


def test_gibbs_and_laplace_agree_with_nuts(hump_data):
    truth, x, y = hump_data
    fits = {s: BADLM(nlag=20, k=12, sampler=s).fit(x, y, samples=400, warmup=400, chains=2, seed=2)
            for s in ("nuts", "gibbs", "laplace")}
    for s, m in fits.items():
        _check_recovery(m, truth)
    ref = fits["nuts"].lag_curve_.reshape(-1, 21)
    sd = ref.std(axis=0)
    # Gibbs targets the same posterior; Laplace is an approximation and is held
    # to a looser standard.
    for s, tol in (("gibbs", 0.5), ("laplace", 1.0)):
        other = fits[s].lag_curve_.reshape(-1, 21)
        assert np.max(np.abs(other.mean(axis=0) - ref.mean(axis=0)) / sd) < tol, s
    assert abs(np.log(fits["gibbs"].draws_["sigma2"].mean() / fits["nuts"].draws_["sigma2"].mean())) < 0.15


@pytest.mark.parametrize("tau_model", ["basis", "monotone", "horseshoe", "constant"])
def test_tau_models_run(hump_data, tau_model):
    truth, x, y = hump_data
    m = BADLM(nlag=20, k=12, tau_model=tau_model).fit(x, y, samples=200, warmup=200, chains=1, seed=3)
    assert np.all(np.isfinite(m.lag_curve_))
    _check_recovery(m, truth, rmse_frac=0.3, min_cover=0.5)


@pytest.mark.parametrize("end_constraint", ["none", "zero_slope"])
def test_end_constraints_run(hump_data, end_constraint):
    truth, x, y = hump_data
    m = BADLM(nlag=20, k=12, end_constraint=end_constraint).fit(x, y, samples=200, warmup=200, chains=1, seed=4)
    _check_recovery(m, truth, rmse_frac=0.3, min_cover=0.5)
    if end_constraint == "zero_slope":
        assert np.allclose(m.lag_curve_[..., -1], 0.0, atol=1e-10)


def test_ar1_errors_recover_phi():
    truth = readme_lag_curve(20)
    x, y = simulate_dlm(truth, n=500, standardize_x=True, ar_eps=0.6, sd_eps=0.01, seed=5)
    m = BADLM(nlag=20, k=12, ar1=True).fit(x, y, samples=300, warmup=300, chains=1, seed=5)
    phi = m.draws_["phi"].ravel()
    assert np.quantile(phi, 0.025) < 0.6 < np.quantile(phi, 0.975)
    _check_recovery(m, truth, rmse_frac=0.3, min_cover=0.6)


def test_poisson_family(rng):
    truth = 10 * readme_lag_curve(20)
    x = rng.normal(size=500)
    eta = 2.0 + lag_matrix(x, 20) @ truth
    y = rng.poisson(np.exp(eta))
    m = BADLM(nlag=20, k=12, family="poisson").fit(x, y, samples=300, warmup=300, chains=1, seed=6)
    _check_recovery(m, truth, rmse_frac=0.3, min_cover=0.6)
    assert abs(m.draws_["intercept"].mean() - 2.0) < 0.2
    pr = m.predict(x[:60], kind="response")
    assert pr["mean"].shape == (40,) and np.all(pr["mean"] > 0)
    assert m.effective_df() is None


def test_covariates_and_offset(rng):
    truth = readme_lag_curve(15)
    x, y = simulate_dlm(truth, n=300, standardize_x=True, ar_eps=0.0, sd_eps=0.01, seed=8)
    z = rng.normal(size=300)
    y = y + 0.05 * z[15:]
    m = BADLM(nlag=15, k=10).fit(x, y, covariates=z, samples=200, warmup=200, chains=1, seed=8)
    gamma = m.draws_["theta"][..., -1].ravel() * m.std_.y_sd
    assert abs(gamma.mean() - 0.05) < 0.01
    _check_recovery(m, truth, rmse_frac=0.3, min_cover=0.6)


def test_invalid_configurations():
    with pytest.raises(ValueError):
        BADLM(nlag=10, sampler="gibbs", family="poisson")
    with pytest.raises(ValueError):
        BADLM(nlag=10, ar1=True, family="poisson")
    with pytest.raises(ValueError):
        BADLM(nlag=10, end_constraint="flat")
    with pytest.raises(ValueError):
        BADLM(nlag=10, k=3)
