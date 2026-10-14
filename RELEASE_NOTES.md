# BioEntropy V1.0

Final desktop companion software for multivariate biomedical data analysis.

## Software

- Normal-reference Mahalanobis deviation uses Ledoit–Wolf full covariance. Group covariance entropy and the separate distance-derived state score are recorded with method identifiers.
- Multi-visit experiments ordinarily export eight main artifacts. Figures 02/04 show aligned nonzero baselines with shared group/point scaling; raw estimates remain in figures 01/03 and workbooks. Single-visit experiments retain nine outputs and their original group comparisons.
- Outlier flags are retained for quality control and do not trigger subject attenuation or exclusion. Run audit records identify inputs, settings, method policies and computational environment.
- The separate forest module computes standardized effects and exploratory cross-model summaries with documented preprocessing.

## Downloads

- **BioEntropy_V1.0_Windows.zip** — extract and run `BioEntropy_V1_0.exe` on Windows x64.
- **BioEntropy_V1.0_Source.zip** — source, public tests, input templates and technical usage documentation.
- **BioEntropy_V1.0_SHA256.txt** — checksums for the archives and executable.

The journal manuscript and Supporting Information are supplied with the article and are not included in this repository or these downloads.

## Verification

The validated calculation and figure behavior are retained. The numerical regression suite passed 162 tests with 7 optional external-data tests skipped. Source and Windows self-tests preserve raw distance and entropy summaries. The tested experiment-106 results match all 10 preceding workbooks, including 142,108 numerical cells; study data are not distributed.
