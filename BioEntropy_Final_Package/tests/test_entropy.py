"""Unit tests for the system-entropy estimator.

compute_group_entropy estimates the differential entropy of a multivariate
Gaussian via a Ledoit-Wolf shrunk covariance:

    H = 0.5 * logdet(Sigma) + 0.5 * k * log(2*pi*e)

We avoid asserting exact theoretical values (shrinkage biases finite samples)
and instead pin invariant properties the science relies on.
"""
import numpy as np
import pytest

from bioentropy_core_legacy import (
    compute_group_entropy,
    fill_nan_with_col_median,
)


def _gaussian(n, k, scale, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, scale, size=(n, k))


def test_returns_finite_for_valid_matrix():
    h = compute_group_entropy(_gaussian(200, 3, 1.0))
    assert np.isfinite(h)


def test_nan_when_fewer_than_two_rows():
    assert np.isnan(compute_group_entropy(np.array([[1.0, 2.0, 3.0]])))


def test_scaling_adds_k_log_c():
    # H(c*X) = H(X) + k*log(c): Ledoit-Wolf shrinkage is scale-equivariant,
    # so scaling every feature by c shifts entropy by exactly k*log(c).
    X = _gaussian(400, 4, 1.0, seed=1)
    k = X.shape[1]
    c = 2.0
    h1 = compute_group_entropy(X)
    h2 = compute_group_entropy(c * X)
    assert h2 - h1 == pytest.approx(k * np.log(c), rel=1e-6, abs=1e-6)


def test_higher_variance_means_higher_entropy():
    low = compute_group_entropy(_gaussian(300, 3, 1.0, seed=2))
    high = compute_group_entropy(_gaussian(300, 3, 5.0, seed=2))
    assert high > low


def test_fill_nan_with_col_median_replaces_only_nans():
    X = np.array([[1.0, np.nan], [3.0, 4.0], [5.0, 6.0]])
    filled = fill_nan_with_col_median(X)
    assert not np.isnan(filled).any()
    # the NaN in column 1 becomes that column's median (of 4 and 6 -> 5)
    assert filled[0, 1] == pytest.approx(5.0)
    # non-NaN entries are untouched
    assert filled[1, 0] == 3.0 and filled[2, 1] == 6.0
