# BioEntropy tests

Run from `BioEntropy_Final_Package/`:

```bash
pip install -r requirements_dev.txt
pytest tests/ -q
```

## What is covered

| File | Scope |
|------|-------|
| `test_hedges_g.py` | Hedges' g formula, J small-sample correction, sign/antisymmetry, edge cases |
| `test_infer_direction.py` | Disease-burden sign convention, incl. the `non-HDL` harmful-vs-beneficial regression guard |
| `test_entropy.py` | System-entropy estimator invariants (scale-equivariance, monotonicity), NaN handling |
| `test_regression_pipeline.py` | End-to-end pipeline output pinned to a committed baseline |

## Regression baseline

`_regression_snapshot.py` runs the real desktop pipeline
(`load_input_bundle → analyze_all → build_markdown_report → generate_output_bundle`)
on the bundled animal-experiment data and hashes the numeric result tables, the
markdown report, and the output-bundle structure into `baseline_snapshot.json`.

`test_regression_pipeline.py` asserts the current output still matches that
baseline. This is the safety net that makes core-module refactors (e.g. dead-code
removal) verifiable. If you make an **intentional** scientific change, regenerate
the baseline:

```bash
python tests/_regression_snapshot.py tests/baseline_snapshot.json
```

and review the diff before committing. The test skips automatically if the raw
data folder is absent (e.g. stripped for a public release).
