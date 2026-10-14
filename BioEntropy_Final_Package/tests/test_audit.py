"""Advisory QC must expose legacy pitfalls without changing scientific results."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import bioentropy_audit as A
import bioentropy_core_generalized as G
import bioentropy_core_legacy as L
import bootstrap_run as B


def test_capture_hashes_the_bytes_actually_loaded(tmp_path):
    path = tmp_path / "SYN-0W.xlsx"
    original = pd.DataFrame({"Group": ["A", "B"], "Sample ID": [1, 2], "f0": [1.5, 2.5]})
    original.to_excel(path, index=False)
    uploads, records = A.capture_input_files([path])
    path.write_bytes(b"changed after capture")
    raw, *_ = L.load_input_bundle(uploads)
    pd.testing.assert_frame_equal(raw[("SYN", 0)], original)
    assert records[0]["sha256"] == hashlib.sha256(uploads[0][1]).hexdigest()


def test_missing_nmd_advisory_preserves_zero_score(tmp_path):
    sample = pd.DataFrame({"Week": [0], "Group": ["Treatment"], "Sample ID": [1],
                           "NMD": [0.0], "NMD_feature_count": [0]})
    before = sample.copy(deep=True)
    result = {"experiments": {"SYN": {"sample_level": sample, "true_normal_available": True,
                                       "true_normal_params": pd.DataFrame({"Feature": ["f0", "f1"]})}}}
    output = tmp_path / "Result"
    output.mkdir()
    (output / "existing.csv").write_text("NMD\n0\n", encoding="utf-8")
    existing = (output / "existing.csv").read_bytes()
    bundle = A.write_audit_bundle(result, output)
    pd.testing.assert_frame_equal(sample, before)
    assert (output / "existing.csv").read_bytes() == existing
    audit = Path(bundle["audit_dir"])
    rows = pd.read_csv(audit / "qc_advisories.csv")
    assert "NMD_NO_OBSERVED_FEATURES" in rows["code"].tolist()
    manifest = json.loads((audit / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["input_capture"] == "not_provided"
    assert manifest["winsor_q"] is None


def test_unit_and_normal_fallback_are_advisory_only():
    params = pd.DataFrame({"Feature": ["glucose (mg/dL)"], "Normal_max": [6.1]})
    before = params.copy(deep=True)
    results = {"experiments": {"MOUSE": {"builtin_reference_range_used": True,
                                          "external_range_params": params,
                                          "true_normal_available": False,
                                          "score_mode": "builtin_reference_range"}}}
    findings = A.collect_advisories(results, {"normal_reference_experiments": ["MOUSE"]})
    codes = {record["code"] for record in findings}
    assert {"BUILTIN_ADULT_RANGE_CONTEXT", "BUILTIN_RANGE_UNIT_REVIEW", "UPLOADED_NORMAL_FALLBACK"} <= codes
    pd.testing.assert_frame_equal(params, before)


def test_audit_serializes_numpy_and_nonfinite_values(tmp_path):
    output = tmp_path / "Result"
    output.mkdir()
    params = pd.DataFrame({"Feature": ["f0"], "value": [np.nan], "count": [np.int64(3)]})
    results = {"experiments": {"SYN": {"preprocess": params}}}
    audit = A.write_audit_bundle(results, output, {"winsor_q": np.float64(0.01)})
    data = json.loads((Path(audit["audit_dir"]) / "run_manifest.json").read_text(encoding="utf-8"))
    assert data["experiments"]["SYN"]["preprocessing"][0]["value"] is None
    assert data["experiments"]["SYN"]["preprocessing"][0]["count"] == 3


def test_audit_failure_keeps_completed_legacy_bundle(tmp_path, monkeypatch):
    output = tmp_path / "Result"
    output.mkdir()
    (output / "existing.csv").write_text("value\n42\n", encoding="utf-8")
    expected = (output / "existing.csv").read_bytes()
    monkeypatch.setattr(G._legacy, "generate_output_bundle", lambda *a, **k: {
        "output_dir": str(output), "zip_path": str(tmp_path / "Result.zip")})
    monkeypatch.setattr(G, "_append_joint_visual_artifacts", lambda *a: None)
    import bioentropy_figure_data as F
    monkeypatch.setattr(F, "export_figure_data_tables", lambda *a: [])

    def fail(*args):
        raise OSError("audit disk unavailable")

    monkeypatch.setattr(A, "write_audit_bundle", fail)
    bundle = G.generate_output_bundle({"experiments": {}, "tables": {}})
    assert "audit disk unavailable" in bundle["audit_error"]
    assert Path(bundle["zip_path"]).exists()
    assert (output / "existing.csv").read_bytes() == expected


def test_cleanup_refuses_result_directory_outside_owned_temp_root(tmp_path):
    output = tmp_path / "user_results" / "Result"
    output.mkdir(parents=True)
    assert A.cleanup_generated_bundle({"output_dir": str(output)}) is False
    assert output.exists()


def test_cleanup_owned_root_only_after_copy(tmp_path, monkeypatch):
    monkeypatch.setattr(A.tempfile, "gettempdir", lambda: str(tmp_path))
    root = tmp_path / "bioentropy_run_test"
    output = root / "Result"
    output.mkdir(parents=True)
    assert A.cleanup_generated_bundle({"output_dir": str(output)}) is True
    assert not root.exists()


def test_working_bootstrap_runtime_is_not_reinstalled(monkeypatch):
    monkeypatch.setattr(B, "ensure_venv", lambda: None)
    monkeypatch.setattr(B, "vpy", lambda: Path("python"))
    monkeypatch.setattr(B.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0))
    commands = []
    monkeypatch.setattr(B, "run", commands.append)
    B.main()
    assert commands == [[Path("python"), B.APP]]


@pytest.mark.parametrize("supports_lock, expected", [(True, B.LOCK), (False, B.REQ)])
def test_bootstrap_selects_python_compatible_requirements(monkeypatch, supports_lock, expected):
    monkeypatch.setattr(B, "ensure_venv", lambda: None)
    monkeypatch.setattr(B, "vpy", lambda: Path("python"))
    probes = iter([1, 0 if supports_lock else 1])
    monkeypatch.setattr(B.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=next(probes)))
    commands = []
    monkeypatch.setattr(B, "run", commands.append)
    B.main()
    assert commands[0] == [Path("python"), "-m", "pip", "install", "-r", expected]
    assert all("--upgrade" not in command for command in commands)


def test_exported_template_can_be_loaded(tmp_path):
    path = tmp_path / "TEMPLATE-0W.xlsx"
    L.build_excel_template(path)
    raw, *_ = L.load_input_bundle([path])
    frame = raw[("TEMPLATE", 0)]
    assert frame["Sample ID"].is_unique
    assert len(frame) == 4
