# BioEntropy

BioEntropy is a small, self-contained desktop tool for the exploratory,
multi-index quantification of *disease deviation* and *system entropy* in
grouped biomedical measurements (e.g. clinical laboratory panels or animal
biochemistry). It computes, per group and per time point:

- **NMD** — a robust Normal-referenced Mahalanobis deviation (distance of a
  group from a healthy/Normal reference cloud);
- **System entropy** — the log-determinant (Gaussian) entropy of the group's
  multivariate distribution, estimated with Ledoit–Wolf shrinkage;
- Absolute (**NRBS**, published-range) or baseline-relative (**BRI**) burden
  scores when no enrolled Normal arm is available;
- **Hedges' g** effect sizes and random-effects pooled estimates for
  treatment-versus-comparator contrasts.

It is intended as a *supporting analytical utility*, not as a confirmatory
statistical package: its correlation and pooled outputs are exploratory and
should be interpreted alongside a pre-specified analysis plan. The exact
definitions, equations and parameter settings are documented in
[`BioEntropy_Final_Package/METHODS.md`](BioEntropy_Final_Package/METHODS.md).

## Installation

Requires Python 3.11+.

```bash
cd BioEntropy_Final_Package
python -m pip install -r requirements_desktop.txt
```

## Running

Launch the desktop application:

```bash
cd BioEntropy_Final_Package
python BioEntropy_Desktop_App.py
```

Load one or more standardized `.xlsx` workbooks (one column per group; a
separate `Normal.xlsx` may be supplied as a healthy reference), assign the
treatment / comparator / reference roles, and run the analysis. Results
(figures and result tables) are written to a timestamped output folder.

## Tests

```bash
cd BioEntropy_Final_Package
python -m pip install -r requirements_dev.txt
MPLBACKEND=Agg python -m pytest tests/
```

The suite includes unit tests for the numerical primitives and a regression
test that pins the full pipeline output (tables and every figure's pixels) to a
committed baseline.

## Repository layout

- `BioEntropy_Final_Package/` — application source, methods reference, and tests.

## Data and licensing

This project is released under the MIT License (see [`LICENSE`](LICENSE)).

The research measurements used during development are non-public and are not
distributed with this repository; supply your own standardized workbooks (one
column per group, one row per subject) to run the tool. With no datasets
present, the regression test self-skips while the unit tests run in full.
