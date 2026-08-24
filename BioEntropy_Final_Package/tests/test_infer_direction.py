"""Unit tests for infer_direction (disease-burden sign convention).

Convention:
    +1  -> a higher value usually indicates WORSE status (harmful marker)
    -1  -> a higher value usually indicates BETTER status (beneficial marker)

These tests pin the documented bug fix: 'non-HDL' / '非高密度脂蛋白' must be
treated as harmful (+1) even though it contains the beneficial 'HDL' substring.
"""
import pytest

from bioentropy_core_generalized import infer_direction


@pytest.mark.parametrize(
    "name",
    ["LDL-C", "ldl", "总胆固醇", "甘油三酯", "TG", "空腹血糖", "glucose",
     "HbA1c", "糖化血红蛋白", "CRP", "ALT", "AST", "GGT"],
)
def test_harmful_markers_positive(name):
    assert infer_direction(name) == 1.0


@pytest.mark.parametrize(
    "name",
    ["HDL-C", "高密度脂蛋白", "高密度脂蛋白胆固醇", "ApoA1", "载脂蛋白A1",
     "脂联素", "adiponectin"],
)
def test_beneficial_markers_negative(name):
    assert infer_direction(name) == -1.0


@pytest.mark.parametrize("name", ["非高密度脂蛋白", "non-HDL", "nonHDL", "Non-HDL-C"])
def test_non_hdl_is_harmful_not_beneficial(name):
    # Regression guard: the harmful 'non-HDL' must win over the 'HDL' substring.
    assert infer_direction(name) == 1.0


def test_unknown_marker_defaults_to_harmful():
    # Unknown markers default to +1 (conservative: treat as disease-burden).
    assert infer_direction("某未知指标XYZ") == 1.0


def test_case_and_whitespace_insensitive():
    assert infer_direction("  hdl-c  ") == infer_direction("HDL-C")
