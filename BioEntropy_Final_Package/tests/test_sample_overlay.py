"""Tests for the per-subject point overlay on single-timepoint group bars.

The bars are unchanged; the overlay reads the render-pass sample-level frame from
display state and maps each metric to its per-subject values (NMD/NRBS/BRI
directly; Entropy via the reference-state contribution 0.5·ln(2πe·NMD²)).
"""
import numpy as np
import pandas as pd
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bioentropy_display_state as _dstate
from bioentropy_plots import _per_sample_points_for_metric, _overlay_trajectory_points


@pytest.fixture
def sample_level():
    return pd.DataFrame({
        "Experiment": ["E1"] * 6,
        "Group": ["Normal", "Normal", "Model", "Model", "Drug", "Drug"],
        "NMD": [1.0, 2.0, 10.0, 12.0, 5.0, 6.0],
        "BRI": [0.1, 0.2, 0.9, 1.0, 0.5, 0.4],
    })


def _with_state(sl, fn):
    prev = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    _dstate.CURRENT_SAMPLE_LEVEL = sl
    try:
        return fn()
    finally:
        _dstate.CURRENT_SAMPLE_LEVEL = prev


def test_nmd_maps_to_raw_per_sample(sample_level):
    pts = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "NMD_mean", ["Normal", "Model", "Drug"]))
    assert set(pts) == {"Normal", "Model", "Drug"}
    assert sorted(pts["Model"].tolist()) == [10.0, 12.0]


def test_entropy_uses_reference_state_contribution(sample_level):
    pts = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "Entropy", ["Normal", "Model"]))
    # 0.5·ln(2πe·NMD²) is monotone increasing in NMD, so Model dots sit above Normal
    assert float(np.mean(pts["Model"])) > float(np.mean(pts["Normal"]))
    expected = 0.5 * np.log(2 * np.pi * np.e * 10.0 ** 2)
    assert abs(float(pts["Model"][0]) - expected) < 1e-9


def test_no_state_means_no_points(sample_level):
    prev = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    _dstate.CURRENT_SAMPLE_LEVEL = None
    try:
        assert _per_sample_points_for_metric("E1", "NMD_mean", ["Normal"]) == {}
    finally:
        _dstate.CURRENT_SAMPLE_LEVEL = prev


def test_unknown_metric_returns_empty(sample_level):
    pts = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "NHPS_mean", ["Normal"]))
    assert pts == {}


def _weekly_sample_level():
    rows = []
    for wk in (0, 6, 12):
        for grp, base in [("A", 2.0), ("B", 8.0)]:
            for i in range(5):
                rows.append({"Experiment": "E1", "Week": wk, "Group": grp,
                             "NMD": float(base + i * 0.2 + wk * 0.01)})
    return pd.DataFrame(rows)


def test_trajectory_overlay_plots_all_weeks():
    sl = _weekly_sample_level()
    fig, ax = plt.subplots()
    cmap = {"A": "#123456", "B": "#654321"}
    vals = _with_state(sl, lambda: _overlay_trajectory_points(
        ax, "E1", "NMD_mean", [0, 6, 12], ["A", "B"], cmap))
    plt.close(fig)
    # 2 groups × 3 weeks × 5 subjects = 30 plotted points
    assert len(vals) == 30
    assert min(vals) >= 2.0 and max(vals) <= 9.0


def test_trajectory_overlay_normalized_is_0_100():
    sl = _weekly_sample_level()
    fig, ax = plt.subplots()
    cmap = {"A": "#123456", "B": "#654321"}
    vals = _with_state(sl, lambda: _overlay_trajectory_points(
        ax, "E1", "NMD_mean", [0, 6, 12], ["A", "B"], cmap, normalize="display"))
    plt.close(fig)
    assert vals and min(vals) >= 0.0 and max(vals) <= 100.0


def test_trajectory_overlay_no_state_is_empty():
    fig, ax = plt.subplots()
    prev = getattr(_dstate, "CURRENT_SAMPLE_LEVEL", None)
    _dstate.CURRENT_SAMPLE_LEVEL = None
    try:
        assert _overlay_trajectory_points(ax, "E1", "NMD_mean", [0], ["A"], {"A": "#000"}) == []
    finally:
        _dstate.CURRENT_SAMPLE_LEVEL = prev
        plt.close(fig)


def test_display_normalization_anchors_0_100(sample_level):
    """display normalization anchors Normal→0 and Model→100 (with gamma)."""
    pts = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "NMD_mean", ["Normal", "Model", "Drug"], normalize="display"))
    allv = np.concatenate([pts[g] for g in pts])
    assert allv.min() >= 0.0 and allv.max() <= 100.0
    # Normal subjects sit near 0, Model subjects near 100
    assert float(np.mean(pts["Normal"])) < float(np.mean(pts["Drug"])) < float(np.mean(pts["Model"]))


def test_linear_normalization_differs_from_display(sample_level):
    lin = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "NMD_mean", ["Normal", "Model", "Drug"], normalize="linear"))
    disp = _with_state(sample_level, lambda: _per_sample_points_for_metric(
        "E1", "NMD_mean", ["Normal", "Model", "Drug"], normalize="display"))
    # gamma>1 pushes mid values down, so Drug's display mean is below its linear mean
    assert float(np.mean(disp["Drug"])) < float(np.mean(lin["Drug"]))
