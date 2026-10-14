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

[Download BioEntropy V1.0](https://github.com/MrElss/BioEntropy/releases/tag/V1.0)
for the Windows application, source archive, usage documentation and checksums.

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
separately in the technical methods document.

Tables and figures are exported with a sibling audit directory containing input
hashes, software and method identifiers, environment information and advisory QC.
Outlier flags do not trigger subject-level attenuation or exclusion. Other defined
preprocessing and missing-data rules still apply.

Multi-visit experiments output eight main artifacts, with no duplicate figures 09/10.
Single-visit experiments retain their nine artifacts and original normalized
figures 02/04. Raw estimates remain in figures 01/03 and result workbooks.

Longitudinal figures 02 and 04 align each group's baseline before shared
0–100 normalization and gamma=3 display. Each group subtracts its own baseline
and adds the arithmetic mean of eligible group baselines. The aligned group
estimates determine one common range. Baseline points therefore coincide at a
nonzero display value; decreasing or unchanged trajectories start at 100.
When there is no decrease and the range is nonconstant, the lower anchor is
extended downward by the observed span. A constant range uses linear 50.
Missing 0W uses the labelled earliest valid visit; at least two valid visits are
required. Burden points share their group's offset and the curve's anchors.
Entropy has no individual dots. Raw estimates, differences and display anchors
remain available. This index is not a raw metric or percent change, and zero
does not mean zero distance or entropy. It does not add paired inference or a
baseline-adjusted treatment effect.

Normalized figures use the same group-summary reference anchors for the group
curves/bars and subject-level disease-deviation points. Clipping and the cubic
display transform still mean that a plotted group value need not equal the mean
of its plotted points. Entropy is displayed at group level; log-distance subject
points are not overlaid as individual entropy observations. Raw metrics and
statistical results are unchanged by these display rules.

The final GitHub release is **BioEntropy V1.0**. The software is a companion
tool for journal research. The manuscript and journal Supporting Information
are maintained separately and are not distributed with this repository or release.

Each study must document its sample sizes, selected markers, reference groups,
settings and contrasts. Composite pooling across different marker panels needs
a comparability argument. Distance-derived scores are mathematically related
to NMD and are not independent validation.

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
