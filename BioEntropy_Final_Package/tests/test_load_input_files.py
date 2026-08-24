"""Guard for the load_input_files / run_from_paths entry point.

load_input_files unpacks the result of load_input_bundle. load_input_bundle
returns a 7-tuple, so the unpacking must stay in sync (a previous mismatch
silently broke the public run_from_paths entry point with a ValueError).
"""
import glob
from pathlib import Path

import pytest

import bioentropy_core_generalized as G

DATA = Path(__file__).resolve().parent.parent.parent / "data" / "动物实验"


def _animal_files():
    return sorted(glob.glob(str(DATA / "**" / "*.xlsx"), recursive=True))


def test_load_input_bundle_returns_seven_tuple():
    # Pin the arity that load_input_files depends on.
    files = _animal_files()
    if not files:
        pytest.skip("bundled animal data not present")
    result = G.load_input_bundle(files)
    assert len(result) == 7


@pytest.mark.skipif(not _animal_files(), reason="bundled animal data not present")
def test_load_input_files_unpacks_cleanly():
    raw, overview = G.load_input_files(_animal_files())
    assert isinstance(raw, dict) and len(raw) > 0


@pytest.mark.skipif(not _animal_files(), reason="bundled animal data not present")
def test_run_from_paths_end_to_end(tmp_path, monkeypatch):
    # run_from_paths writes an output bundle to the working directory.
    monkeypatch.chdir(tmp_path)
    results = G.run_from_paths(_animal_files())
    assert results.get("tables")
    assert "bundle" in results
