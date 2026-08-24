"""Tests for upfront experiment-file validation (_validate_week_dataframe)."""
import glob
from pathlib import Path

import pandas as pd
import pytest

import bioentropy_core_legacy as L

DATA = Path(__file__).resolve().parent.parent.parent / "data" / "动物实验"


def _good():
    return pd.DataFrame({
        "Group": ["A", "A", "B", "B"],
        "Sample ID": [1, 2, 3, 4],
        "glucose": [5.1, 5.4, 9.2, 8.8],
        "ldl": [2.1, 2.3, 3.9, 4.1],
    })


def test_valid_frame_passes():
    L._validate_week_dataframe("ok.xlsx", _good())  # must not raise


def test_duplicate_sample_id_raises():
    df = _good()
    df.loc[3, "Sample ID"] = 1  # duplicate of row 0
    with pytest.raises(ValueError, match="重复的 Sample ID"):
        L._validate_week_dataframe("dup.xlsx", df)


def test_all_text_feature_column_raises():
    df = _good()
    df["glucose"] = ["low", "high", "mid", "n/a"]
    with pytest.raises(ValueError, match="没有任何可识别的数值"):
        L._validate_week_dataframe("text.xlsx", df)


def test_empty_feature_column_is_tolerated():
    df = _good()
    df["spare"] = [None, None, None, None]
    L._validate_week_dataframe("empty_col.xlsx", df)  # must not raise


def test_non_numeric_sample_ids_are_tolerated():
    # Clinical files use patient codes; upstream coerces them to NaN. Validation
    # must NOT reject this (regression guard for an over-strict earlier check).
    df = _good()
    df["Sample ID"] = [float("nan")] * 4
    L._validate_week_dataframe("clinical_like.xlsx", df)  # must not raise


def test_no_feature_columns_raises():
    df = _good()[["Group", "Sample ID"]]
    with pytest.raises(ValueError, match="没有任何指标列"):
        L._validate_week_dataframe("nofeat.xlsx", df)


@pytest.mark.skipif(not glob.glob(str(DATA / "**" / "*.xlsx"), recursive=True),
                    reason="bundled data absent")
def test_bundled_data_passes_validation():
    # the shipped real files must all pass (no false positives)
    files = sorted(glob.glob(str(DATA / "**" / "*.xlsx"), recursive=True))
    raw, *_ = L.load_input_bundle(files)
    assert len(raw) > 0
