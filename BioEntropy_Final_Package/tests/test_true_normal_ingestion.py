"""Contract tests for ingesting an uploaded true-Normal sample file.

Covers the three things that let a bare ``Normal.xlsx`` (an ID column + marker
columns, in conventional US units) drive the NMD path for a clinical study
whose data are in SI units:

  * a prefix-less Normal/健康 workbook is recognised as a GLOBAL Normal reference
    and attached to every experiment in the batch;
  * an ID-only layout (no Group column) is parsed as one Normal cloud;
  * conventional (mg/dL) glucose/lipid values are auto-converted to SI (mmol/L)
    so the healthy reference lands on the study's scale.

If the ingestion rules change intentionally, update this test alongside the code.
"""
import pandas as pd
import pytest

import bioentropy_core_legacy as L


def test_bare_normal_filename_is_global_reference():
    assert L.parse_input_file_kind("Normal.xlsx") == ("normal_like", None, None)
    assert L.parse_input_file_kind("正常.xlsx")[0] in ("normal_like", "week")  # 正常 not in alias set → week is fine
    # a hyphen-prefixed normal file stays experiment-scoped
    assert L.parse_input_file_kind("105-Normal.xlsx") == ("normal_like", "105", None)
    # an ordinary week file is untouched
    assert L.parse_input_file_kind("105-12W.xlsx") == ("week", "105", 12)


def test_id_only_normal_sample_parses_as_one_cloud():
    df = pd.DataFrame({
        "健康人ID": [1010410301, 1010410402, 1010410901, 1010411101],
        "血糖": [90.3, 80.6, 94.6, 96.6],
        "总胆固醇": [161.2, 148.8, 151.9, 127.1],
        "甘油三酯": [134.5, 76.1, 110.6, 132.7],
    })
    parsed = L._parse_true_normal_sample_df(df, "Normal.xlsx")
    assert list(parsed.columns[:2]) == ["Group", "Sample ID"]
    assert set(parsed["Group"].unique()) == {"Normal"}
    # no marker column was consumed as a group / ID
    assert "血糖" in parsed.columns and "总胆固醇" in parsed.columns


def test_units_auto_convert_mgdl_to_mmoll():
    # study cholesterol in mmol/L (~4.4); normal in mg/dL (~155) → factor 0.02586
    study = pd.Series([4.2, 4.4, 4.6, 5.0])
    normal = pd.Series([150.0, 155.0, 160.0, 145.0])
    conv, factor, note = L._harmonize_true_normal_units("总胆固醇 (mmol/L)-105-ADEFF", normal, study)
    assert factor == pytest.approx(0.02586, rel=1e-3)
    assert note == "auto_converted_mgdl_to_mmoll"
    assert conv.median() == pytest.approx(normal.median() * 0.02586, rel=1e-3)


def test_units_left_alone_when_already_matched():
    # both in mmol/L; healthy < diseased but same unit → no conversion
    study = pd.Series([2.4, 2.6, 2.8])
    normal = pd.Series([2.0, 2.2, 2.3])
    _conv, factor, note = L._harmonize_true_normal_units("低密度脂蛋白胆固醇 (mmol/L)-105-ADEFF", normal, study)
    assert factor == 1.0
    assert note == ""


def test_bare_glucose_gets_semantic_code_without_colliding():
    assert L.feature_semantic_code("血糖 (mmol/L)-105-ADLB") == "glu_blood"
    assert L.feature_semantic_code("血糖") == "glu_blood"
    # qualified glucose and HbA1c keep their own codes
    assert L.feature_semantic_code("空腹血浆血糖 (mmol/L)-105-ADEFF") == "glu_fasting"
    assert L.feature_semantic_code("糖化血红蛋白（HbA1c）") == "hba1c"
