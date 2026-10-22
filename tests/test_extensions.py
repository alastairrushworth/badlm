"""Multi-exposure curves, observation weights and modifier (tensor-product) surfaces."""
import jax.numpy as jnp
import numpy as np
import pytest

from badlm import BADLM, lag_matrix
from badlm.collapsed import component_precision
from badlm.design import LagComponent, SuffStats, build_modifier, build_structure


@pytest.mark.parametrize("ec", ["zero", "none"])
def test_kronecker_sum_logdet_matches_eigenvalues(ec, rng):
    st = build_structure(20, 10, end_constraint=ec)
    comp = LagComponent(st, build_modifier(rng.normal(size=300), K=5))
    tau = rng.normal(size=st.K_tau)
    Q, ld = component_precision(comp, jnp.asarray(tau), 2.3)
    ev = np.linalg.eigvalsh(np.asarray(Q))
    assert np.isclose(float(ld), np.sum(np.log(ev[ev > 1e-8])))
    assert np.sum(ev <= 1e-8) == (1 if ec == "none" else 0)


def test_weighted_sufficient_statistics(rng):
    X, y, w = rng.normal(size=(40, 3)), rng.normal(size=40), rng.uniform(0.5, 2.0, 40)
    a = SuffStats.from_data(X, y, weights=w)
    b = SuffStats.from_data(X * np.sqrt(w)[:, None], y * np.sqrt(w))
    for name in ("S_a", "S_b", "S_c", "s_a", "s_b", "s_c"):
        assert np.allclose(getattr(a, name), getattr(b, name))
    assert np.isclose(a.yy_b, b.yy_b)


@pytest.fixture(scope="module")
def two_exposure_data():
    rng = np.random.default_rng(7)
    n, lag = 2500, np.arange(21.0)
    x1, x2 = rng.normal(size=n), rng.normal(size=n)
    b1 = 0.5 * np.exp(-lag / 3.0)
    b2 = -0.4 * np.exp(-0.5 * ((lag - 8) / 2.5) ** 2)
    b1[-1] = b2[-1] = 0.0
    mu = lag_matrix(x1, 20) @ b1 + lag_matrix(x2, 20) @ b2
    return x1, x2, mu + rng.normal(size=mu.size), b1, b2, mu


def test_two_exposures_recovered_with_laplace(two_exposure_data):
    x1, x2, y, b1, b2, mu = two_exposure_data
    m = BADLM(nlag=[20, 20], k=[14, 14], sampler="laplace").fit(np.c_[x1, x2], y, samples=300, chains=1)
    assert np.sqrt(np.mean((m.lag_curve_summary(exposure=0).mean - b1) ** 2)) < 0.03
    assert np.sqrt(np.mean((m.lag_curve_summary(exposure=1).mean - b2) ** 2)) < 0.03
    assert m.draws_["theta"].shape[-1] == 2 * m.structs_[0].K_r + 1
    pr = m.predict(np.c_[x1, x2][-100:])
    assert pr["mean"].shape == (80,)
    assert np.corrcoef(pr["mean"], mu[-80:])[0, 1] > 0.98
    assert "lag_curve_1" in m.diagnostics_


def test_modifier_surface_recovers_amplitude(two_exposure_data):
    x1, _, _, b1, _, _ = two_exposure_data
    rng = np.random.default_rng(3)
    z = np.abs(rng.normal(size=x1.size)) + 0.2
    amp = 1.0 + 1.5 * (z - z.mean())
    mu = (lag_matrix(x1, 20) * amp[20:, None]) @ b1
    y = mu + rng.normal(size=mu.size)
    m = BADLM(nlag=20, k=14, modifier_k=5, sampler="laplace").fit(x1, y, modifier=z, samples=300, chains=1)
    zs = np.quantile(z, [0.2, 0.8])
    surf = m.lag_surface(zs)
    for i, zz in enumerate(zs):
        truth = b1 * (1.0 + 1.5 * (zz - z.mean()))
        assert np.sqrt(np.mean((surf["mean"][i] - truth) ** 2)) < 0.05
    # the surface at the average state is close to the average-amplitude curve
    assert np.sqrt(np.mean((m.summary_.mean - b1) ** 2)) < 0.05
    pr = m.predict(x1[-60:], modifier=z[-60:])
    assert pr["mean"].shape == (40,)


def test_time_modifier_and_weights_run(two_exposure_data):
    x1, _, y, _, _, _ = two_exposure_data
    mt = BADLM(nlag=20, k=14, modifier_k=4, sampler="laplace").fit(x1, y, modifier="time", samples=100, chains=1)
    assert mt.lag_surface([0.0, 1.0])["median"].shape == (2, 21)
    mw = BADLM(nlag=20, k=14, sampler="laplace").fit(x1, y, weights=np.full(y.size, 2.0), samples=100, chains=1)
    assert np.isfinite(mw.summary_.median).all()
    with pytest.raises(ValueError):
        BADLM(nlag=20, k=14, modifier_k=4).fit(x1, y)
    with pytest.raises(ValueError):
        BADLM(nlag=[20, 20], sampler="gibbs")


def test_log_lag_scale_and_temper(two_exposure_data):
    x1, _, y, b1, _, _ = two_exposure_data
    from badlm.priors import Priors

    m = BADLM(nlag=20, k=12, lag_scale="log", sampler="laplace", priors=Priors(temper=4.0)).fit(
        x1, y, samples=100, chains=1)
    assert m.struct.B.shape == (21, 12)
    assert np.isfinite(m.summary_.median).all()
    with pytest.raises(ValueError):
        BADLM(nlag=20, lag_scale="sqrt")


def test_walk_forward_predictions_are_out_of_sample(two_exposure_data):
    from badlm.backtest import performance, trade, walk_forward

    x1, _, y, b1, _, mu = two_exposure_data
    time = np.arange(x1.size)
    model = BADLM(nlag=20, k=10)
    wf = walk_forward(model, x1, y, time, start=1200, step=400)
    assert wf.mean.size == x1.size - 1 - 1200          # rows with target time after the first refit
    assert np.all(wf.time > 1200)
    # predictions use only earlier data yet track the part of the mean due to x1
    mu1 = lag_matrix(x1, 20) @ b1
    assert np.corrcoef(wf.mean, mu1[wf.row - 20])[0, 1] > 0.9
    tr = trade(wf.mean, wf.var, np.zeros(wf.mean.size), wf.time, rule="gp", speed=0.5)
    perf = performance(tr["position"], tr["trade"], wf.y, wf.time, costs=(0.0, 0.01))
    assert perf["sharpe_gross"] > 0 and perf["sharpe_net"][0.01] < perf["sharpe_gross"]


def test_exposure_surface_recovers_nonlinear_response():
    rng = np.random.default_rng(5)
    n, lag = 4000, np.arange(21.0)
    x = rng.normal(size=n)
    b = 0.5 * np.exp(-lag / 4.0)
    b[-1] = 0.0
    L = lag_matrix(x, 20)
    mu = (L * (1 - 0.3 * L**2)) @ b                      # large moves respond less
    y = mu + rng.normal(size=mu.size)
    m = BADLM(nlag=20, k=12, modifier_k=5, sampler="laplace").fit(x, y, modifier="exposure", samples=100, chains=1)
    s = m.lag_surface([-1.0, 0.0, 1.0])
    for i, xv in enumerate(s["z"]):
        assert abs(s["median"][i, 0] - b[0] * (1 - 0.3 * xv**2)) < 0.1
    pr = m.predict(x[-300:])
    assert np.corrcoef(pr["mean"], mu[-280:])[0, 1] > 0.9


def test_ar1_carry_over_improves_walk_forward_log_score():
    from badlm.backtest import log_score, walk_forward
    from badlm.simulate import ar1

    rng = np.random.default_rng(11)
    n, lag = 3000, np.arange(21.0)
    x = rng.normal(size=n)
    b = 0.3 * np.exp(-lag / 4.0)
    b[-1] = 0.0
    mu = lag_matrix(x, 20) @ b
    y = mu + ar1(mu.size, 0.7, 1.0, rng)                      # strongly autocorrelated errors
    time = np.arange(n)
    plain = walk_forward(BADLM(nlag=20, k=10), x, y, time, start=1500, step=500)
    with_ar = walk_forward(BADLM(nlag=20, k=10, ar1=True), x, y, time, start=1500, step=500)
    assert with_ar.fit_info[-1]["phi"] > 0.5
    assert log_score(with_ar.y, with_ar.mean, with_ar.var) > log_score(plain.y, plain.mean, plain.var) + 0.1


def test_blocked_cv_and_ensemble(two_exposure_data):
    from badlm.backtest import blocked_cv, ensemble, subperiod_sharpe

    x1, _, y, b1, _, _ = two_exposure_data
    time = np.arange(x1.size)
    cv = blocked_cv(BADLM(nlag=20, k=10), x1, y, time, n_blocks=4, embargo=20)
    assert cv.mean.size == x1.size - 20
    mu1 = lag_matrix(x1, 20) @ b1
    assert np.corrcoef(cv.mean, mu1[cv.row - 20])[0, 1] > 0.9
    ens = ensemble([cv, cv])
    assert np.allclose(ens.mean, cv.mean) and np.allclose(ens.var, cv.var)
    sp = subperiod_sharpe(np.sign(cv.mean) * cv.y, cv.time, breaks=[1000, 2000])
    assert len(sp) == 3 and all(s[2] > 0 for s in sp)



def test_panel_model_recovers_loadings_and_curves():
    from badlm.panel import PanelFitter, build_panel_design, panel_walk_forward

    rng = np.random.default_rng(3)
    N, T, p = 8, 1500, 15
    lag = np.arange(p + 1.0)
    g1 = 0.4 * np.exp(-lag / 3.0)
    g1[-1] = 0.0
    g2 = 0.3 * np.exp(-0.5 * ((lag - 5) / 2.0) ** 2)
    g2[-1] = 0.0
    f = rng.normal(size=T)
    u1, u2 = 1 + 0.5 * rng.normal(size=N), 1 + 0.3 * rng.normal(size=N)
    ys, fs, xs, assets, times, mus = [], [], [], [], [], []
    for i in range(N):
        x = rng.normal(size=T)
        mu = u1[i] * (lag_matrix(f, p) @ g1) + u2[i] * (lag_matrix(x, p) @ g2)
        ys.append(np.r_[np.zeros(p), mu] + rng.normal(size=T))
        fs.append(f); xs.append(x); assets.append(np.full(T, i)); times.append(np.arange(T)); mus.append(mu)
    y, F, X, A, Tm = map(np.concatenate, (ys, fs, xs, assets, times))
    d = build_panel_design([F, X], y, Tm, A, nlag=p, k=10)
    st = d.stats(np.ones(d.y.size, bool))
    fit = PanelFitter(d, st)
    res = fit.fit(st)
    c1, c2 = fit.curves(res)
    # the scale is shared between a curve and its loadings; the products are identified
    prod1 = res["loading"][:, 0][:, None] * c1[None, :]
    prod2 = res["loading"][:, 1][:, None] * c2[None, :]
    assert np.sqrt(np.mean((prod1 - u1[:, None] * g1[None, :]) ** 2)) < 0.05
    assert np.sqrt(np.mean((prod2 - u2[:, None] * g2[None, :]) ** 2)) < 0.05
    assert np.corrcoef(res["loading"][:, 0], u1)[0, 1] > 0.95
    mu_true = np.concatenate(mus)
    wf = panel_walk_forward(d, start=900, step=300)
    wfp = panel_walk_forward(d, start=900, step=300, pool=True)
    sel = np.isin(d.row, wf["row"])
    assert np.corrcoef(wf["mean"], mu_true[sel])[0, 1] > np.corrcoef(wfp["mean"], mu_true[sel])[0, 1] + 0.03
    # modifier surfaces inside the panel run
    z = rng.normal(size=y.size)
    dm = build_panel_design([F, X], y, Tm, A, nlag=p, k=10, modifier=z, modifier_k=4)
    stm = dm.stats(np.ones(dm.y.size, bool))
    resm = PanelFitter(dm, stm).fit(stm)
    assert np.isfinite(resm["mean"]).all() and len(PanelFitter(dm, stm).curves(resm)[0]) == 5
