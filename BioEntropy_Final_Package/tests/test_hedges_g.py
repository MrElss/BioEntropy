"""Unit tests for the Hedges' g effect-size calculation.

Validates the documented formulas (see PROJECT_NOTES.md / METHODS.md):

    pooled SD  = sqrt(((n1-1)*SD1^2 + (n2-1)*SD2^2) / (n1+n2-2))
    Cohen's d  = (mean1 - mean2) / pooled SD
    J          = 1 - 3 / (4*(n1+n2-2) - 1)
    Hedges' g  = J * Cohen's d

The implementation now lives in ``gvalue_recompute`` (the old
``bioentropy_hedges_forest`` module was removed); ``hedges_g`` returns a
``GResult`` dataclass.
"""
import math

import numpy as np
import pytest

from gvalue_recompute import hedges_g


def _manual_g(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n1, n2 = len(x), len(y)
    sd1, sd2 = x.std(ddof=1), y.std(ddof=1)
    df = n1 + n2 - 2
    pooled = math.sqrt(((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / df)
    d = (x.mean() - y.mean()) / pooled
    j = 1.0 - 3.0 / (4.0 * df - 1.0)
    return pooled, d, j * d


def test_matches_manual_formula():
    x = [5.0, 6.0, 7.0, 8.0, 9.0]
    y = [1.0, 2.0, 3.0, 4.0, 5.0]
    pooled, d, g = _manual_g(x, y)
    res = hedges_g(x, y)
    assert res is not None
    assert res.pooled_sd == pytest.approx(pooled, rel=1e-12)
    assert res.d == pytest.approx(d, rel=1e-12)
    assert res.g == pytest.approx(g, rel=1e-12)


def test_se_and_ci_are_consistent():
    res = hedges_g([5, 6, 7, 8, 9], [1, 2, 3, 4, 5])
    assert res.se > 0
    assert res.ci_low < res.g < res.ci_high
    # CI is the documented normal approximation g +/- z*SE.
    z = (res.ci_high - res.ci_low) / (2 * res.se)
    assert z == pytest.approx(1.959963984540054, rel=1e-9)


def test_j_correction_shrinks_toward_zero():
    # |g| must be strictly smaller than |d| (small-sample bias correction).
    res = hedges_g([10, 12, 14, 16], [1, 2, 3, 4])
    assert abs(res.g) < abs(res.d)
    df = res.s1.n + res.s2.n - 2
    j = 1.0 - 3.0 / (4.0 * df - 1.0)
    assert res.g == pytest.approx(j * res.d, rel=1e-12)


def test_sign_follows_mean_difference():
    higher = hedges_g([10, 11, 12, 13], [1, 2, 3, 4])
    lower = hedges_g([1, 2, 3, 4], [10, 11, 12, 13])
    assert higher.g > 0
    assert lower.g < 0
    # antisymmetry of swapping the two samples
    assert higher.g == pytest.approx(-lower.g, rel=1e-12)


def test_identical_samples_zero_effect():
    res = hedges_g([4, 5, 6, 7], [4, 5, 6, 7])
    assert res.g == pytest.approx(0.0, abs=1e-12)


def test_too_few_samples_returns_none():
    assert hedges_g([1.0], [1.0, 2.0, 3.0]) is None


def test_zero_variance_returns_none():
    # pooled SD == 0 -> undefined effect size
    assert hedges_g([5, 5, 5, 5], [5, 5, 5, 5]) is None


def test_handles_nan_entries():
    # NaNs should be dropped, not propagate into the result.
    res = hedges_g([5, 6, 7, 8, float("nan")], [1, 2, 3, 4])
    assert res is not None
    assert math.isfinite(res.g)
    # the NaN row is dropped, leaving n=4 in the first group
    assert res.s1.n == 4
