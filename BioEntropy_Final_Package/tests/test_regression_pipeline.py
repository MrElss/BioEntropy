"""End-to-end regression test.

Runs the real desktop pipeline (load_input_bundle -> analyze_all ->
build_markdown_report -> generate_output_bundle) on the bundled animal-experiment
data and asserts the numeric/structural output is byte-identical to a committed
baseline. This is the safety net that lets the core modules be refactored
(e.g. dead-code removal) with confidence.

If the raw data folder is absent (e.g. stripped for a public release), the test
skips rather than fails.
"""
import json
import os
from pathlib import Path

import pytest

import _regression_snapshot as snap

BASELINE = Path(__file__).resolve().parent / "baseline_snapshot.json"


@pytest.mark.skipif(not snap._animal_files(), reason="bundled animal data not present")
def test_pipeline_output_matches_baseline():
    assert BASELINE.exists(), "baseline_snapshot.json missing; regenerate with _regression_snapshot.py"

    # This test pins exact table/report/figure BYTES, which are environment
    # specific: figure pixels depend on the installed fonts (CJK labels), on
    # whether plotly is present, and on the matplotlib backend, not just on
    # library versions. It is therefore a same-machine refactoring guard, not a
    # portable check. Skip it on CI (where fonts/plotly differ); the rest of the
    # suite is environment-robust. Force it locally with BIOENTROPY_BYTE_REGRESSION=1.
    if os.environ.get("CI") and not os.environ.get("BIOENTROPY_BYTE_REGRESSION"):
        pytest.skip("byte/pixel regression is environment-specific (fonts, plotly, backend); skipped on CI")

    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    baseline_versions = baseline.pop(snap.VERSION_KEY, None)

    # Secondary guard for non-CI runs on a different numeric/plotting stack.
    current_versions = snap.env_versions()
    if baseline_versions is not None and baseline_versions != current_versions:
        pytest.skip(
            "library versions differ from the baseline-recording environment, so "
            "the byte/pixel regression does not apply here.\n"
            f"  baseline: {baseline_versions}\n  current:  {current_versions}\n"
            "Regenerate with `python tests/_regression_snapshot.py` to re-pin."
        )

    current = snap.compute_snapshot()
    # Compare per-key for a readable failure listing which table/report drifted.
    drifted = sorted(k for k in set(baseline) | set(current)
                     if baseline.get(k) != current.get(k))
    assert not drifted, f"output changed for: {drifted}"
