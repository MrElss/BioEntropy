"""Unit tests for the added statistical methodology:
Benjamini-Hochberg FDR and the entropy reliability diagnostics.
"""
import numpy as np
import pytest

from bioentropy_entropy import (
    benjamini_hochberg,
    compute_group_entropy,
    compute_group_entropy_with_diagnostics,
)


# ---------------------------------------------------------------- BH FDR ----
def test_bh_matches_hand_computation():
    p = [0.001, 0.01, 0.02, 0.5]
    q = benjamini_hochberg(p)
    # q_i = p_(i) * m / rank, then monotone from the top
    assert q[0] == pytest.approx(0.004, rel=1e-9)   # 0.001*4/1
    assert q[1] == pytest.approx(0.02, rel=1e-9)    # 0.01*4/2
    assert q[2] == pytest.approx(0.02 + 2/300, rel=1e-6)  # 0.02*4/3 = 0.02667
    assert q[3] == pytest.approx(0.5, rel=1e-9)     # 0.5*4/4


def test_bh_q_at_least_p_and_bounded():
    p = np.array([0.04, 0.03, 0.2, 0.7, 0.001])
    q = benjamini_hochberg(p)
    assert np.all(q >= p - 1e-12)          # q-values never below raw p
    assert np.all((q >= 0) & (q <= 1))     # bounded to [0, 1]


def test_bh_is_monotone_in_p_order():
    p = np.array([0.001, 0.004, 0.008, 0.02, 0.6])  # already ascending
    q = benjamini_hochberg(p)
    assert np.all(np.diff(q) >= -1e-12)    # non-decreasing with p


def test_bh_preserves_input_order_and_nan():
    q = benjamini_hochberg([0.5, np.nan, 0.001])
    assert np.isnan(q[1])                  # NaN stays NaN, in place
    assert q[2] < q[0]                     # smaller p -> smaller q
    # NaN entries are excluded from the test count m (m=2 here)
    assert q[2] == pytest.approx(0.002, rel=1e-9)  # 0.001*2/1


def test_bh_empty_and_all_nan():
    assert benjamini_hochberg([]).shape == (0,)
    assert np.all(np.isnan(benjamini_hochberg([np.nan, np.nan])))


# ----------------------------------------------- entropy diagnostics ----
def _gauss(n, p, seed=0):
    return np.random.default_rng(seed).normal(size=(n, p))


def test_diagnostics_value_matches_plain_entropy():
    X = _gauss(50, 4, seed=1)
    h_plain = compute_group_entropy(X)
    h_diag, _ = compute_group_entropy_with_diagnostics(X)
    assert h_diag == pytest.approx(h_plain, rel=1e-12)


def test_diagnostics_report_shape_and_ratio():
    X = _gauss(30, 6, seed=2)
    _, d = compute_group_entropy_with_diagnostics(X)
    assert d["entropy_n"] == 30 and d["entropy_p"] == 6
    assert d["entropy_n_over_p"] == pytest.approx(5.0)
    assert 0.0 <= d["entropy_ledoitwolf_shrinkage"] <= 1.0


@pytest.mark.parametrize("n,p,expected", [
    (60, 6, "reliable"),   # n/p = 10
    (18, 6, "limited"),    # n/p = 3
    (8, 16, "low"),        # n < p  (the animal-data regime)
])
def test_reliability_bins(n, p, expected):
    _, d = compute_group_entropy_with_diagnostics(_gauss(n, p, seed=3))
    assert d["entropy_cov_reliability"] == expected
