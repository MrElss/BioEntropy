"""Contract tests for the published-reference-range → NRBS absolute anchor.

When a study has no enrolled Normal/healthy arm and no uploaded Normal-range
file, clinical-lab markers are anchored to built-in published adult reference
ranges so an *absolute* disease-deviation score (NRBS) replaces the purely
baseline-relative BRI. These tests pin:

  * the marker→range matching (longest, most specific token wins; direction is
    encoded by which bound is open), and
  * the figure metric preference NMD > NRBS > BRI.

If the reference table or preference order is changed intentionally, update this
test together with the code so the contract stays honest.
"""
import math

import numpy as np
import pandas as pd
import pytest

import bioentropy_core_legacy as L
from bioentropy_plots import _preferred_burden_kind


def _clinical_raw():
    """Minimal 2-week, 2-group clinical-style panel with standard lab columns."""
    cols = [
        "空腹血浆血糖 (mmol/L)-105-ADEFF",
        "低密度脂蛋白胆固醇 (mmol/L)-105-ADEFF",
        "非高密度脂蛋白胆固醇 (mmol/L)-105-ADEFF",
        "高密度脂蛋白胆固醇 (mmol/L)-105-ADEFF",
        "丙氨酸氨基转移酶 (U/L)-105-ADEFF",
    ]
    raw = {}
    rng = np.random.default_rng(0)
    for w in (0, 12):
        n = 24
        data = {"Group": ["HTD1801"] * (n // 2) + ["安慰剂"] * (n // 2),
                "Sample ID": list(range(n))}
        for c in cols:
            data[c] = rng.normal(5.0, 1.0, n)
        raw[("105", w)] = pd.DataFrame(data)
    return raw


def test_builtin_ranges_match_by_longest_token():
    raw = _clinical_raw()
    rdf = L.build_default_reference_range_df(raw, "105")
    by_feat = {r["Feature"]: r for _, r in rdf.iterrows()}

    # Higher-worse labs: only an upper bound is set (lower is open / NaN).
    fpg = by_feat["空腹血浆血糖 (mmol/L)-105-ADEFF"]
    assert math.isclose(fpg["Normal_max"], 6.1) and pd.isna(fpg["Normal_min"])

    # Non-HDL must bind to the non-HDL limit (4.1), NOT the HDL limit — the
    # longer / more specific token wins even though 「高密度脂蛋白胆固醇」 is a substring.
    nonhdl = by_feat["非高密度脂蛋白胆固醇 (mmol/L)-105-ADEFF"]
    assert math.isclose(nonhdl["Normal_max"], 4.1) and pd.isna(nonhdl["Normal_min"])

    # HDL-C is higher-is-better: only a lower bound is set.
    hdl = by_feat["高密度脂蛋白胆固醇 (mmol/L)-105-ADEFF"]
    assert math.isclose(hdl["Normal_min"], 1.0) and pd.isna(hdl["Normal_max"])


def test_builtin_ranges_empty_for_unmatched_panel():
    raw = {("A", 0): pd.DataFrame({"Group": ["m", "c"], "Sample ID": [0, 1],
                                   "weird_marker_x": [1.0, 2.0],
                                   "weird_marker_y": [3.0, 4.0]})}
    assert L.build_default_reference_range_df(raw, "A").empty


def test_analyze_experiment_uses_builtin_ranges_and_computes_nrbs():
    raw = _clinical_raw()
    res = L.analyze_experiment(raw, "105", treatment_group="HTD1801", comparator_group="安慰剂")
    assert res["score_mode"] == "builtin_reference_range"
    assert res["builtin_reference_range_used"] is True
    assert res["true_normal_available"] is False
    assert pd.to_numeric(res["summary"]["NRBS_mean"], errors="coerce").notna().any()


def test_preferred_burden_kind_prefers_nrbs_over_bri():
    # No true Normal, but NRBS available → NRBS beats BRI.
    summary = pd.DataFrame({"BRI_mean": [0.1, -0.2], "NRBS_mean": [0.3, 0.1]})
    assert _preferred_burden_kind({"summary": summary, "true_normal_available": False}) == "NRBS"


def test_preferred_burden_kind_prefers_nmd_over_nrbs():
    summary = pd.DataFrame({"BRI_mean": [0.1], "NRBS_mean": [0.3], "NMD_mean": [4.0]})
    assert _preferred_burden_kind({"summary": summary, "true_normal_available": True}) == "NMD"


def test_preferred_burden_kind_falls_back_to_bri():
    summary = pd.DataFrame({"BRI_mean": [0.1, -0.2], "NRBS_mean": [np.nan, np.nan]})
    assert _preferred_burden_kind({"summary": summary, "true_normal_available": False}) == "BRI"


@pytest.mark.parametrize("token,feature,expected_max", [
    ("糖化血红蛋白", "糖化血红蛋白 (%)-106-ADLB", 6.0),
    ("甘油三酯", "甘油三酯 (mmol/L)-106-ADEFF", 1.7),
    ("总胆固醇", "总胆固醇 (mmol/L)-106-ADEFF", 5.2),
])
def test_builtin_ranges_cover_common_labs(token, feature, expected_max):
    raw = {("106", 0): pd.DataFrame({"Group": ["m", "c"], "Sample ID": [0, 1], feature: [5.0, 6.0]}),
           ("106", 12): pd.DataFrame({"Group": ["m", "c"], "Sample ID": [0, 1], feature: [5.0, 6.0]})}
    rdf = L.build_default_reference_range_df(raw, "106")
    assert not rdf.empty
    assert math.isclose(float(rdf.iloc[0]["Normal_max"]), expected_max)
