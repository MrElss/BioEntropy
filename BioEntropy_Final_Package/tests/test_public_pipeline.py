"""Always-run public regression: numerical values, workbook export, figures and audit isolation."""
import json
import zipfile
from pathlib import Path

import pytest

import bioentropy_core_generalized as G
from bioentropy_audit import capture_input_files
import _public_regression as R
from bioentropy_constants import SOFTWARE_VERSION


@pytest.fixture(scope="module", params=R.CASES)
def public_case(request, tmp_path_factory):
    return request.param, R.run_case(tmp_path_factory.mktemp("public_inputs"), request.param)


def _compare(actual, expected, location):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), location
        for key in expected:
            _compare(actual[key], expected[key], f"{location}/{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), location
        for index, (a, e) in enumerate(zip(actual, expected)):
            _compare(a, e, f"{location}/{index}")
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, rel=1e-9, abs=1e-10), location
    else:
        assert actual == expected, location


def test_public_numeric_baseline(public_case):
    case, (results, *_) = public_case
    baseline = json.loads(R.BASELINE.read_text(encoding="utf-8"))
    assert baseline["source_commit"] == R.SOURCE_COMMIT
    assert baseline["calculation_version"] == SOFTWARE_VERSION
    _compare(R.table_snapshot(results), baseline["cases"][case], case)


def test_public_export_and_audit_do_not_change_results(public_case):
    case, (results, files, _raw, _ranges, normals, roles) = public_case
    before = R.table_snapshot(results)
    _, inputs = capture_input_files(files)
    bundle = G.generate_output_bundle(
        results, package_name="PublicRegression",
        audit_context={"inputs": inputs, "role_map": roles, "winsor_q": G.WINSOR_Q_DEFAULT,
                       "normal_reference_experiments": list(normals)},
    )
    output, audit = Path(bundle["output_dir"]), Path(bundle["audit_dir"])
    assert audit.parent == output.parent and audit != output
    assert R.table_snapshot(results) == before
    layout = json.loads(Path(__file__).with_name("public_artifact_layout.json").read_text(encoding="utf-8"))
    assert sorted(R.artifact_snapshot(output)) == layout[case]
    assert len(bundle["figure_paths"]) > 0
    assert len(bundle["figure_data_files"]) > 0
    manifest = json.loads((audit / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["inputs"] == inputs
    assert manifest["winsor_q"] == G.WINSOR_Q_DEFAULT
    assert manifest["method_policy"]["primary_nmd"] == "ledoitwolf_full_covariance"
    assert manifest["method_policy"]["primary_nmd_z_cap"] is None
    assert manifest["experiments"][case]["roles"]["treatment"] == "Treatment"
    with zipfile.ZipFile(bundle["zip_path"]) as archive:
        assert all(name.startswith("PublicRegression/") for name in archive.namelist())
        assert not any("run_manifest" in name or "qc_advisories" in name for name in archive.namelist())
