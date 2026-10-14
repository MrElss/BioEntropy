# BioEntropy V1.0 tests

Run `python -m pip install -r requirements_dev.txt` and `python -m pytest tests/ -q`
from the application directory. `ruff check .` runs the configured static checks.

The public suite includes analytical NMD checks, covariance-entropy invariants,
Hedges' g formulas, input validation, reference ranges, quality-control behavior,
fixed synthetic numeric regression, workbook and figure exports, output layout,
and audit isolation. `test_calculation_methods.py` checks the primary and diagonal
sensitivity routes against an independent synthetic reference fixture.

`public_numeric_baseline.json` defines the current synthetic numeric contract.
`reference_metric_baseline.json` is an independent synthetic comparison fixture.
Do not regenerate numeric fixtures to conceal an unexpected change. Review the
method and the resulting differences before accepting any intentional revision.

Some input-specific checks require external research files and skip when those
files are absent. Research data and private regression snapshots are not distributed.
The source and Windows executable also support `--self-test <report.json>` using
synthetic inputs to exercise loading, analysis, table/figure export and run records.
