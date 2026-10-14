"""Read-only provenance and advisory QC, written beside the legacy result bundle.

This module does not fit estimators, resample, change input data, or add columns
to scientific result tables. QC findings describe the existing computation.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata as metadata
import json
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib
from matplotlib import font_manager
import numpy as np
import pandas as pd

from bioentropy_constants import SOFTWARE_VERSION, RNG_SEED, NMD_PRIMARY_METHOD, ENTROPY_SUBJECT_POLICY, NMD_Z_CAP

PKG = Path(__file__).resolve().parent
LIBRARIES = ("numpy", "pandas", "scipy", "scikit-learn", "matplotlib",
             "openpyxl", "xlsxwriter", "plotly")
QC_COLUMNS = ("code", "level", "experiment", "week", "group", "sample_id", "detail")


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_input_files(paths: Iterable[str | Path]) -> tuple[list[tuple[str, bytes]], list[dict]]:
    """Hash exactly the workbook bytes that will be handed to the existing loader."""
    uploads, records = [], []
    for item in paths:
        path = Path(item)
        content = path.read_bytes()
        uploads.append((path.name, content))
        records.append({"name": path.name, "size_bytes": len(content),
                        "sha256": hashlib.sha256(content).hexdigest()})
    return uploads, records


def _json_value(value: Any):
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if value is None or value is pd.NA:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (str, Path)):
        return str(value)
    raise TypeError(f"Unsupported audit value: {type(value).__name__}")


def _records(frame: Any) -> list[dict]:
    return _json_value(frame.to_dict("records")) if isinstance(frame, pd.DataFrame) else []


def _source_provenance() -> dict:
    files = {}
    for path in sorted(PKG.glob("*.py")):
        files[path.name] = _file_sha256(path)
    build_info = PKG / "BUILD_INFO.json"
    info = json.loads(build_info.read_text(encoding="utf-8")) if build_info.exists() else {}
    commit, dirty = None, None
    if not getattr(sys, "frozen", False):
        flags = subprocess.CREATE_NO_WINDOW if sys.platform == "win32" else 0
        try:
            sha = subprocess.run(["git", "-C", str(PKG), "rev-parse", "HEAD"],
                                 capture_output=True, text=True, timeout=3, creationflags=flags)
            status = subprocess.run(["git", "-C", str(PKG), "status", "--porcelain"],
                                    capture_output=True, text=True, timeout=3, creationflags=flags)
            if sha.returncode == 0:
                commit = sha.stdout.strip()
            if status.returncode == 0:
                dirty = bool(status.stdout.strip())
        except (OSError, subprocess.TimeoutExpired):
            pass
    executable_hash = _file_sha256(Path(sys.executable)) if getattr(sys, "frozen", False) else None
    return {"git_commit": commit, "git_dirty": dirty, "build_info": info,
            "executable_sha256": executable_hash,
            "source_sha256": files, "frozen_executable": bool(getattr(sys, "frozen", False))}


def collect_advisories(results: dict, context: dict | None = None) -> list[dict]:
    """Collect observations only; all current scores and inclusion rules stay intact."""
    findings = []
    context = context or {}

    def add(code, experiment, detail, row=None, level="review"):
        row = row if row is not None else {}
        findings.append({"code": code, "level": level, "experiment": str(experiment),
                         "week": row.get("Week"), "group": row.get("Group"),
                         "sample_id": row.get("Sample ID"), "detail": detail})

    for experiment, result in results.get("experiments", {}).items():
        sample = result.get("sample_level", pd.DataFrame())
        summary = result.get("summary", pd.DataFrame())
        if "NMD_feature_count" in sample and result.get("true_normal_available"):
            count = pd.to_numeric(sample["NMD_feature_count"], errors="coerce")
            params = result.get("true_normal_params", pd.DataFrame())
            p = len(params) if isinstance(params, pd.DataFrame) else 0
            for _, row in sample.loc[count == 0].iterrows():
                add("NMD_NO_OBSERVED_FEATURES", experiment,
                    "全部 NMD 指标缺失；现有算法以 Normal 均值填补，可能得到 NMD=0。该样本未被排除，不能据此解释为健康。", row)
            if p:
                for _, row in sample.loc[(count > 0) & (count / p < 0.7)].iterrows():
                    add("NMD_LOW_FEATURE_COVERAGE", experiment,
                        f"NMD 有效指标不足参考指标的 70%（参考指标数 {p}）；缺失值填补和原分数保持不变。", row)
            add("REFERENCE_P_VALUES_REVIEW", experiment,
                "NMD 的 chi-square/F p 值受有限参考、收缩协方差和预处理影响；请作为参考性输出核查，旁路 QC 未重算 p 值。")
        if "Entropy_method" in summary and summary["Entropy_method"].eq("normal_reference_state_entropy").any():
            add("STATE_ENTROPY_DISTANCE_COUPLING", experiment,
                "主 Entropy 使用 0.5*log(mean(NMD^2))+常数，与 NMD 存在结构性耦合；相关性不能作为独立证据。协方差熵保留在 Entropy_robust_covariance。")
        if "Entropy_cov_reliability" in summary:
            for _, row in summary.loc[summary["Entropy_cov_reliability"].eq("low")].iterrows():
                add("COVARIANCE_ENTROPY_LOW_N_OVER_P", experiment,
                    "协方差熵的 n/p 小于 2，已有可靠性标签为 low；数值及其区间保持不变。", row)
        if result.get("builtin_reference_range_used"):
            add("BUILTIN_ADULT_RANGE_CONTEXT", experiment,
                "采用内置成人参考范围；请核对物种、人群和单位。现有算法未根据物种选择参考范围。")
            for record in _records(result.get("external_range_params")):
                feature = str(record.get("Feature", ""))
                if re.search(r"mg\s*/\s*dl", feature, flags=re.I) and any(
                    token in feature.lower() for token in
                    ("glucose", "血糖", "cholesterol", "胆固醇", "triglyceride", "甘油三酯", "hdl", "ldl", "tg", "tc")
                ):
                    add("BUILTIN_RANGE_UNIT_REVIEW", experiment,
                        f"{feature} 声明 mg/dL，而匹配的内置参考值采用 mmol/L 数值；现有 NRBS 路径未自动换算该单位。")
        if str(experiment) in context.get("normal_reference_experiments", []) and not result.get("true_normal_available"):
            add("UPLOADED_NORMAL_FALLBACK", experiment,
                f"上传的 Normal 参考未形成可用 NMD；实际模式为 {result.get('score_mode')}，请核对指标匹配及样本数。")
        for record in _records(result.get("true_normal_params")):
            if record.get("Unit_factor_applied") not in (None, 1.0):
                add("NORMAL_UNIT_CONVERSION_RECORDED", experiment,
                    f"{record.get('Feature')}: Normal 单位转换因子 {record.get('Unit_factor_applied')}；{record.get('Unit_note')}", level="info")
    for issue in context.get("export_issues", []):
        add("OPTIONAL_EXPORT_FAILED", "", str(issue))
    return _json_value(findings)


def write_audit_bundle(results: dict, output_dir: Path, context: dict | None = None) -> dict:
    """Write a sibling directory; legacy result files and ZIP contents are untouched."""
    context = context or {}
    output_dir = Path(output_dir)
    audit_dir = output_dir.with_name(output_dir.name + "_audit")
    audit_dir.mkdir(parents=True, exist_ok=True)
    versions = {}
    for library in LIBRARIES:
        try:
            versions[library] = metadata.version(library)
        except metadata.PackageNotFoundError:
            versions[library] = None
    experiments = {}
    for experiment, result in results.get("experiments", {}).items():
        experiments[str(experiment)] = {
            "weeks": result.get("weeks"), "anchor_week": result.get("anchor_week"),
            "features": result.get("features"), "score_mode": result.get("score_mode"),
            "roles": {key: result.get(key + "_group") for key in ("treatment", "comparator", "reference")},
            "preprocessing": _records(result.get("preprocess")),
            "normal_reference": _records(result.get("true_normal_params")),
            "reference_ranges": _records(result.get("external_range_params")),
        }
    artifacts = {path.relative_to(output_dir).as_posix(): _file_sha256(path)
                 for path in sorted(output_dir.rglob("*")) if path.is_file()}
    font_path = Path(font_manager.findfont(font_manager.FontProperties(family=matplotlib.rcParams["font.family"])))
    findings = collect_advisories(results, context)
    manifest = {
        "schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "software_version": SOFTWARE_VERSION, "code": _source_provenance(),
        "method_policy": {"primary_nmd": NMD_PRIMARY_METHOD,
                          "primary_nmd_z_cap": None,
                          "diagonal_sensitivity_z_cap": NMD_Z_CAP,
                          "entropy_subject_policy": ENTROPY_SUBJECT_POLICY},
        "environment": {"python": sys.version, "platform": platform.platform(), "libraries": versions,
                        "matplotlib_backend": matplotlib.get_backend(),
                        "font_family": list(matplotlib.rcParams["font.family"]),
                        "resolved_font": {"name": font_path.name, "sha256": _file_sha256(font_path)}},
        "input_capture": "provided" if "inputs" in context else "not_provided",
        "inputs": context.get("inputs", []), "requested_roles": context.get("role_map"),
        "winsor_q": context.get("winsor_q"), "rng_seed": RNG_SEED,
        "display_options": context.get("display_options"), "experiments": experiments,
        "legacy_artifact_sha256": artifacts, "advisory_count": len(findings),
    }
    (audit_dir / "run_manifest.json").write_text(
        json.dumps(_json_value(manifest), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    with (audit_dir / "qc_advisories.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=QC_COLUMNS)
        writer.writeheader()
        writer.writerows(findings)
    dictionary = PKG / "METRIC_DICTIONARY.md"
    if dictionary.exists():
        shutil.copyfile(dictionary, audit_dir / dictionary.name)
    (audit_dir / "README.md").write_text(
        "# 运行记录与质量提示\n\n"
        "此目录位于分析结果目录旁，不写入原结果 ZIP。所有提示仅供核查，未修改、重算或排除数据。\n\n"
        "- `run_manifest.json`：实际输入哈希、代码与环境、分组角色、预处理和参考参数、输出文件哈希。\n"
        "- `qc_advisories.csv`：现有算法的缺失值、参考范围和指标解释提示。\n"
        "- `METRIC_DICTIONARY.md`：指标定义和解释边界。\n\n"
        "文件哈希用于识别文件；Excel/ZIP 时间戳及绘图环境会影响文件字节。"
        "只有在相同环境下通过回归比较，才能确认旧结果保持一致。\n", encoding="utf-8")
    return {"audit_dir": str(audit_dir), "audit_advisory_count": len(findings)}


def cleanup_generated_bundle(bundle: dict) -> bool:
    """Remove only a verified BioEntropy-owned temporary root after successful copies."""
    output = Path(bundle["output_dir"]).resolve()
    root = output.parent
    if root.parent != Path(tempfile.gettempdir()).resolve() or not root.name.startswith("bioentropy_run_"):
        return False
    if not output.is_dir():
        return False
    shutil.rmtree(root)
    return True
