# BioEntropy V1.0

BioEntropy is a desktop application for exploratory analysis of grouped biomedical
measurements. It supports animal experiments and clinical data collected at multiple
time points, with tables, figures and traceable run records.

## Measures

- **NMD** measures each subject's deviation from a Normal reference using a
  Ledoit–Wolf full covariance matrix. A capped diagonal distance is exported as
  a sensitivity analysis.
- **Covariance entropy** describes the dispersion of each group's multivariate
  measurements. The separate **distance-derived state score** summarizes the
  second moment of NMD. The `Entropy_method` field identifies which measure is
  selected for the main `Entropy` output.
- **NRBS** uses supplied or built-in reference ranges. **BRI** describes burden
  relative to the selected study baseline when a Normal reference is unavailable.
- The separate forest module calculates **Hedges' g** and exploratory pooled
  effects. Its preprocessing differs from the main workflow. Its per-animal
  entropy-labelled score is a log-distance transformation, not group covariance
  entropy or independent evidence of efficacy.

See [Methods](BioEntropy_Final_Package/METHODS.md), the
[metric dictionary](BioEntropy_Final_Package/METRIC_DICTIONARY.md), and the
[analysis guide](BioEntropy_Final_Package/ANALYSIS_GUIDE.md).

## Windows application

[Download BioEntropy V1.0](https://github.com/MrElss/BioEntropy/releases/tag/v1.0.1)
for the Windows application, source archive, supplementary methods and checksums.

Unpack the Windows distribution and run `BioEntropy_V1_0.exe`.
The interface displays **BioEntropy V1.0**. Analysis results are saved in a dated
`BioEntropy_Results_V1_0_*` directory. Settings are stored in the user's home
directory as `.bioentropy_v1_0_desktop.json`.

## Run from source

Python 3.11 or newer is required. From `BioEntropy_Final_Package`, run
`Start_Desktop.cmd` on Windows, or install the dependencies and start the program:

```bash
python -m pip install -r requirements_desktop.txt
python BioEntropy_Desktop_App.py
```

For the fixed Python 3.12 analysis stack, use `requirements_lock.txt` in a fresh
environment. The Windows build dependencies are listed in
`requirements_build_lock.txt`. Reproduction also depends on the input data,
settings, Python version and installed fonts.

## Input and output

The main workflow accepts `.xlsx` files with columns
`Group | Sample ID | feature1 | feature2 ...`, one row per subject. Numeric subject
identifiers must be unique within each file. Name longitudinal files
`Experiment-0W.xlsx`, `Experiment-6W.xlsx`, and so on. Normal reference samples and
reference ranges use the separately supplied templates. Assign group roles before
analysis and verify the species and units of any reference ranges.

The forest module's individual-endpoint inputs use one column per group; its
multivariate inputs use one row per animal. These input formats are documented
separately in the supplementary methods.

Tables and figures are exported with a sibling audit directory containing input
hashes, software and method identifiers, environment information and advisory QC.
Outlier flags do not trigger subject-level attenuation or exclusion. Other defined
preprocessing and missing-data rules still apply.

Longitudinal outputs include separate baseline-change figures (09 and 10) for
disease deviation and Entropy. Each group subtracts its own value at 0W, or its
earliest valid visit when 0W is unavailable. Baseline values therefore align at
zero; these figures describe changes in group summaries, not paired-subject
changes or a baseline-adjusted treatment effect. Groups with fewer than two valid
visits have no change trajectory. Single-timepoint experiments retain group comparisons.

Normalized figures use the same group-summary reference anchors for the group
curves/bars and subject-level disease-deviation points. Clipping and the cubic
display transform still mean that a plotted group value need not equal the mean
of its plotted points. Entropy is displayed at group level; log-distance subject
points are not overlaid as individual entropy observations. Raw metrics and
statistical results are unchanged by these display rules.

The application name remains BioEntropy V1.0. Release `v1.0.1` and Windows file
version `1.0.1.0` identify this display correction.

## Supplementary methods

- [中文 SI](docs/SI_CN.md)
- [English SI](docs/SI_EN.md)
- Word documents: [中文](docs/BioEntropy_V1.0_SI_CN.docx) / [English](docs/BioEntropy_V1.0_SI_EN.docx)

The supplementary methods describe software calculations. Each study must report
its actual sample sizes, selected markers, reference groups, settings and contrasts.
Composite pooling across different marker panels requires a comparability argument
and sensitivity analyses. Distance-derived scores are mathematically related to
NMD and should not be interpreted as independent validation.

## Validation and build

```bash
python -m pip install -r requirements_dev.txt
python -m ruff check .
python -m pytest tests/ -q
python BioEntropy_Desktop_App.py --self-test selftest.json
```

The suite checks analytical formulas and fixed synthetic numeric results, along
with workbook/figure export and audit behavior. Optional tests requiring external
research data skip when those data are absent. Build the Windows executable with
`Rebuild_EXE.cmd` or `python bootstrap_build.py`.

## Data and license

Research input files, study outputs, local environments and validation logs are
excluded from the public source distribution. The included templates and public
test fixtures contain example or synthetic data.

BioEntropy is distributed under the [MIT License](LICENSE).
