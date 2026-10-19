import numpy as np
import pytest

from badlm.design import (SuffStats, ar1_whitening_matrix, bspline_basis, build_design,
                          build_structure, lag_matrix)


def lag_matrix_reference(rain, p, start_at_zero=True):
    """Transcription of the original R implementation."""
    rows = []
    for n in range(len(rain) - p):
        seg = rain[n:n + p + (1 if start_at_zero else 0)]
        rows.append(seg[::-1])
    return np.array(rows)


def test_lag_matrix_matches_reference(rng):
    x = rng.normal(size=40)
    for p in (1, 3, 10):
        assert np.allclose(lag_matrix(x, p), lag_matrix_reference(x, p))
        assert np.allclose(lag_matrix(x, p, start_at_zero=False), lag_matrix_reference(x, p, False))


def test_lag_matrix_layout_and_validation():
    m = lag_matrix(np.arange(1, 7), 2)
    assert m.shape == (4, 3)
    assert np.array_equal(m[0], [3, 2, 1])
    assert np.array_equal(m[-1], [6, 5, 4])
    with pytest.raises(ValueError):
        lag_matrix(np.arange(6), 6)
    with pytest.raises(ValueError):
        lag_matrix(np.arange(6), 0)


def test_bspline_basis_partition_of_unity():
    B = bspline_basis(np.arange(31.0), nseg=9, deg=3, domain=(0, 30))
    assert B.shape == (31, 12)
    assert np.all(B >= 0)
    assert np.allclose(B.sum(axis=1), 1.0)
    dB = bspline_basis(np.arange(31.0), nseg=9, deg=3, domain=(0, 30), derivative=1)
    assert np.allclose(dB.sum(axis=1), 0.0, atol=1e-10)


@pytest.mark.parametrize("ec,m", [("none", 0), ("zero", 1), ("zero_slope", 2)])
def test_structure_constraints(ec, m):
    st = build_structure(20, 10, end_constraint=ec)
    assert st.m == m and st.K_r == 10 - m
    assert np.allclose(st.Z.T @ st.Z, np.eye(st.K_r))
    if m >= 1:
        assert np.allclose(st.B_r[-1], 0.0, atol=1e-12)
    if m == 2:
        dB = bspline_basis([20.0], nseg=7, deg=3, domain=(0, 20), derivative=1)[0]
        assert np.allclose(dB @ st.Z, 0.0, atol=1e-12)


@pytest.mark.parametrize("ec", ["none", "zero"])
def test_generalised_determinant_closed_form(ec, rng):
    st = build_structure(20, 10, end_constraint=ec)
    for _ in range(5):
        tau = rng.normal(5, 2, st.K_tau)
        ev = np.linalg.eigvalsh(st.lam_matrix(tau))
        ld = np.log(ev[ev > 1e-9 * ev.max()]).sum()
        assert abs(ld - st.logdet_prior(tau)) < 1e-8
    if ec == "none":
        assert np.isclose(st.logdet_const, np.log(10))   # matrix-tree theorem


def test_coefficients_from_increments_respect_constraints(rng):
    st = build_structure(12, 8, end_constraint="zero")
    d = rng.normal(size=st.K - 1)
    b = st.coefficients_from_increments(d)
    assert np.isclose(st.B[-1] @ b, 0.0)
    assert np.allclose(np.diff(b), d)


def test_ar1_sufficient_statistics_match_whitening(rng):
    X = rng.normal(size=(25, 4))
    y = rng.normal(size=25)
    st = SuffStats.from_data(X, y)
    for phi in (-0.4, 0.0, 0.7):
        L = ar1_whitening_matrix(25, phi)
        S, s, yy = st.at(phi)
        assert np.allclose(S, (L @ X).T @ (L @ X))
        assert np.allclose(s, (L @ X).T @ (L @ y))
        assert np.isclose(yy, (L @ y) @ (L @ y))


def test_build_design_response_handling(rng):
    st = build_structure(5, 5)
    x = rng.normal(size=50)
    y = rng.normal(size=50)
    d_full, _ = build_design(x, y, st)
    d_trim, _ = build_design(x, y[5:], st)
    assert d_full.N == 45 and np.allclose(d_full.y, d_trim.y)
    assert d_full.X.shape == (45, st.K_r + 1)
    with pytest.raises(ValueError):
        build_design(x, y[:30], st)
    with pytest.raises(ValueError):
        build_design(np.r_[x, np.nan], np.r_[y, 1.0], st)


def test_build_design_groups_keep_lags_within_group(rng):
    st = build_structure(3, 5)
    x = rng.normal(size=20)
    y = rng.normal(size=20)
    groups = np.repeat([0, 1], 10)
    d, _ = build_design(x, y, st, groups=groups, standardize=False)
    single = [build_design(x[g == groups], y[g == groups], st, standardize=False)[0] for g in (0, 1)]
    assert d.N == 14
    assert np.allclose(d.L, np.vstack([s.L for s in single]))
    assert np.allclose(d.y, np.concatenate([s.y for s in single]))
    with pytest.raises(ValueError):
        build_design(x, y, st, groups=np.r_[np.zeros(17), 1, 1, 1])
