"""Regression snapshot: runs the real GUI pipeline on the bundled animal AND
clinical experiment data and serialises numeric result tables, the markdown
report, and the output-bundle/explanation-report path to a stable, hashable
structure. Used to prove refactors do not change scientific output.
"""
from __future__ import annotations
import glob, hashlib, json, sys, tempfile
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parent.parent
DATA_ROOT = PKG.parent / "data"
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

# Each dataset is (label, glob pattern). Labels namespace the hash keys so a
# drift report points at the dataset that changed.
DATASETS = {
    "animal": str(DATA_ROOT / "动物实验" / "**" / "*.xlsx"),
    "clinical": str(DATA_ROOT / "临床实验105和106" / "**" / "*.xlsx"),
}


def _files(label: str):
    return sorted(glob.glob(DATASETS[label], recursive=True))


def _animal_files():  # kept for callers/tests that reference it
    return _files("animal")


def available_datasets():
    return [label for label in DATASETS if _files(label)]


def _sha(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def _snapshot_one(label: str) -> dict:
    import bioentropy_core_generalized as G
    files = _files(label)
    raw, range_refs, normal_refs, *_ = G.load_input_bundle(files)
    results = G.analyze_all(
        raw, role_map=None, winsor_q=G.WINSOR_Q_DEFAULT,
        external_range_map=range_refs, true_normal_map=normal_refs,
    )
    snap: dict[str, str] = {}
    p = f"{label}::"
    # (a) numeric result tables
    for name, df in sorted(results.get("tables", {}).items()):
        if not isinstance(df, pd.DataFrame):
            continue
        num = df.select_dtypes(include=[np.number]).round(6)
        snap[p + "table::" + str(name)] = _sha(num.to_csv(index=False))
    # (b) markdown report (exercises the build_markdown_report wrapper chain)
    snap[p + "report::markdown"] = _sha(G.build_markdown_report(results))
    # (c) output bundle structure + explanation-report content
    #     (exercises generate_output_bundle and the explanation-report writer)
    import os
    with tempfile.TemporaryDirectory() as td:
        cwd = os.getcwd(); os.chdir(td)
        try:
            bundle = G.generate_output_bundle(results, package_name="RegressionRun")
            sig = {
                "figure_keys": sorted((bundle.get("figure_paths") or {}).keys()),
                "table_keys": sorted(results.get("tables", {}).keys()),
            }
            snap[p + "bundle::structure"] = _sha(json.dumps(sig, ensure_ascii=False))
            out_dir = Path(bundle["output_dir"])
            md_path = (out_dir / "analysis_report.md").resolve()
            if md_path.exists():
                snap[p + "report::explanation_md"] = _sha(md_path.read_text(encoding="utf-8"))
            # (d) pixel-level figure verification: every PNG/HTML output is byte
            #     deterministic (only .xlsx embeds a volatile zip timestamp), so we
            #     hash each figure's bytes keyed by its path relative to the bundle.
            for fp in sorted(out_dir.rglob("*")):
                if fp.is_file() and fp.suffix.lower() in (".png", ".html"):
                    rel = fp.relative_to(out_dir).as_posix()
                    snap[p + "figure::" + rel] = hashlib.sha256(fp.read_bytes()).hexdigest()
        finally:
            os.chdir(cwd)
    return snap


def compute_snapshot() -> dict:
    snap: dict[str, str] = {}
    for label in available_datasets():
        snap.update(_snapshot_one(label))
    return snap


# Exact table/report/figure bytes depend on the numeric and plotting libraries.
# The baseline records the versions it was generated with so the regression test
# can skip (rather than spuriously fail) on a different environment, e.g. CI
# installing newer matplotlib than the baseline was recorded with.
VERSION_KEY = "__versions__"
_LIBS = ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]


def env_versions() -> dict:
    import importlib.metadata as md
    out = {}
    for lib in _LIBS:
        try:
            out[lib] = md.version(lib)
        except Exception:
            out[lib] = "absent"
    return out


if __name__ == "__main__":
    snap = compute_snapshot()
    payload = dict(snap)
    payload[VERSION_KEY] = env_versions()
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else PKG / "tests" / "baseline_snapshot.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {len(snap)} hashes (+versions {payload[VERSION_KEY]}) to {out}")
