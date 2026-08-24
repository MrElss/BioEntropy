"""Regression tests for header-less group files in ``read_groups``.

Some exported per-metric workbooks omit the header row (the first row is
already numeric data instead of the group names Control / Model / HTD1801 / PM).
Before the fix the first data row was consumed as a header, the columns were
labelled by their numeric value, and the endpoint was silently dropped from the
within-model pool. ``read_groups`` now detects a header-less 4-column sheet and
maps the columns positionally to the standard arm order, recording the
assumption as a QC flag.
"""
import openpyxl

from gvalue_recompute import read_groups


def _write(path, rows):
    wb = openpyxl.Workbook()
    ws = wb.active
    for r in rows:
        ws.append(list(r))
    wb.save(path)


def test_headerless_four_columns_mapped_positionally(tmp_path):
    """No header row -> columns become Control/Model/HTD1801/PM, nothing lost."""
    f = tmp_path / "地鼠-体重.xlsx"
    _write(f, [
        [154.0, 188.0, 160.0, 180.0],
        [150.0, 190.0, 162.0, 178.0],
        [158.0, 186.0, 158.0, 182.0],
    ])
    groups, flagged = read_groups(str(f), model="地鼠", endpoint="Body weight")
    assert set(groups) == {"Control", "Model", "HTD1801", "PM"}
    # every data row is retained (no row consumed as a phantom header)
    assert all(len(v) == 3 for v in groups.values())
    assert groups["Control"] == [154.0, 150.0, 158.0]
    assert groups["Model"] == [188.0, 190.0, 186.0]
    # the assumption is surfaced, not silent
    assert any(fl.kind == "headerless" for fl in flagged)


def test_header_row_still_parsed_normally(tmp_path):
    """A normal header row is unaffected and raises no headerless flag."""
    f = tmp_path / "地鼠-体重.xlsx"
    _write(f, [
        ["Control", "Model", "HTD1801", "PM"],
        [154.0, 188.0, 160.0, 180.0],
        [150.0, 190.0, 162.0, 178.0],
    ])
    groups, flagged = read_groups(str(f), model="地鼠", endpoint="Body weight")
    assert set(groups) == {"Control", "Model", "HTD1801", "PM"}
    assert groups["Control"] == [154.0, 150.0]
    assert not any(fl.kind == "headerless" for fl in flagged)


def test_headerless_non_four_columns_not_guessed(tmp_path):
    """Positional mapping is only assumed for the 4-arm layout, never wider."""
    f = tmp_path / "weird.xlsx"
    _write(f, [
        [1.0, 2.0, 3.0],
        [1.1, 2.1, 3.1],
    ])
    groups, flagged = read_groups(str(f), model="X", endpoint="Y")
    # no positional assumption for a non-standard column count
    assert not any(fl.kind == "headerless" for fl in flagged)
    assert "Control" not in groups
