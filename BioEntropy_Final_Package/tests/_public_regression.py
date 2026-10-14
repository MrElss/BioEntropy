"""Public, fixed synthetic inputs and snapshots; no private research data needed.

The numeric baseline was recorded before the compatibility-only improvements.
No random inputs are generated: fixed observations keep fixtures reproducible.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parent.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))
from bioentropy_constants import SOFTWARE_VERSION

BASELINE = Path(__file__).with_name("public_numeric_baseline.json")
CASES = ("nmd_longitudinal", "nrbs_single", "bri_single")
SOURCE_COMMIT = "1c54f3c0f848694e97a5c0a0b743508c0c6a98ca"

# Three independent feature patterns, repeated only in the larger Normal cloud.
OBSERVATIONS = (
    (-1.2, -0.3, 0.8), (-0.8, 0.9, -0.4), (-0.4, -1.1, -0.7),
    (-0.1, 0.4, 1.2), (0.2, -0.6, -1.0), (0.5, 1.3, 0.1),
    (0.8, -0.8, 0.5), (1.0, 0.2, -0.3),
)


def write_inputs(root: Path, case: str) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    features = (["glucose (mmol/L)", "HDL-C (mmol/L)", "ALT (U/L)"]
                if case == "nrbs_single" else ["f0", "f1", "f2"])
    paths = []
    for week in ((0, 6) if case == "nmd_longitudinal" else (0,)):
        rows = []
        for gi, group in enumerate(("Treatment", "Comparator")):
            for i, observation in enumerate(OBSERVATIONS):
                shift = 1.5 if week == 0 else (0.4 if gi == 0 else 1.2)
                if case == "nrbs_single":
                    values = (7.0 + gi + observation[0],
                              1.1 - 0.1 * gi + 0.1 * observation[1],
                              35.0 + 10 * gi + 8 * observation[2])
                else:
                    values = tuple(v + shift + gi * 0.1 for v in observation)
                row = {"Group": group, "Sample ID": gi * 8 + i + 1,
                       **dict(zip(features, values))}
                if case == "nmd_longitudinal" and week == 6 and gi == 0:
                    if i == 0:
                        row.update(dict.fromkeys(features, np.nan))
                    elif i == 1:
                        row[features[0]] = np.nan
                rows.append(row)
        path = root / f"{case}-{week}W.xlsx"
        pd.DataFrame(rows).to_excel(path, index=False)
        paths.append(path)
    if case == "nmd_longitudinal":
        rows = []
        for i in range(24):
            values = [v + 0.05 * (i // 8 - 1) for v in OBSERVATIONS[i % 8]]
            rows.append({"Group": "Normal", "Sample ID": 100 + i,
                         **dict(zip(features, values))})
        path = root / f"{case}-Normal.xlsx"
        pd.DataFrame(rows).to_excel(path, index=False)
        paths.append(path)
    return paths


def run_case(root: Path, case: str):
    import bioentropy_core_generalized as G

    files = write_inputs(root / case, case)
    raw, ranges, normals, *_ = G.load_input_bundle(files)
    roles = {case: {"treatment": "Treatment", "comparator": "Comparator"}}
    results = G.analyze_all(raw, role_map=roles,
                            external_range_map=ranges, true_normal_map=normals)
    return results, files, raw, ranges, normals, roles


def scalar(value):
    if value is None or value is pd.NA:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if math.isnan(value):
            return None
        return value if math.isfinite(value) else str(value)
    return str(value)


def table_snapshot(results: dict) -> dict:
    return {
        name: {"columns": [str(c) for c in frame.columns],
               "index": [scalar(i) for i in frame.index],
               "rows": [[scalar(v) for v in row] for row in frame.itertuples(index=False, name=None)]}
        for name, frame in sorted(results["tables"].items())
        if isinstance(frame, pd.DataFrame)
    }


def artifact_snapshot(output_dir: Path) -> dict:
    """Compare all legacy files, including workbook XML, excluding ZIP timestamps."""
    snap = {}
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(output_dir).as_posix()
        if path.suffix == ".xlsx":
            members = {}
            with zipfile.ZipFile(path) as archive:
                for name in sorted(archive.namelist()):
                    content = archive.read(name)
                    if name == "docProps/core.xml":
                        content = re.sub(
                            rb"(<dcterms:(?:created|modified)[^>]*>).*?(</dcterms:(?:created|modified)>)",
                            rb"\1TIMESTAMP\2", content)
                    members[name] = hashlib.sha256(content).hexdigest()
            snap[relative] = members
        else:
            snap[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return snap


def record_numeric_baseline(root: Path, destination: Path) -> None:
    payload = {"source_commit": SOURCE_COMMIT, "calculation_version": SOFTWARE_VERSION, "cases": {}}
    for case in CASES:
        results, *_ = run_case(root, case)
        payload["cases"][case] = table_snapshot(results)
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    args = parser.parse_args()
    record_numeric_baseline(args.inputs, args.output)
