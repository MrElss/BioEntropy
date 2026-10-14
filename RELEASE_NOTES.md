# BioEntropy V1.0

## Downloads

- **BioEntropy_V1.0_Windows.zip**: unpack and run `BioEntropy_V1_0.exe` on Windows x64. Includes input templates, methods and Chinese/English supplementary methods.
- **BioEntropy_V1.0_Source.zip**: source, tests, build scripts, templates and documentation. Python 3.11 or newer is required.
- **BioEntropy_V1.0_SI_CN.docx / BioEntropy_V1.0_SI_EN.docx**: editable supplementary methods; each contains 7 pages and 12 editable equations.
- **BioEntropy_V1.0_SHA256.txt**: SHA-256 checksums for the four download files and the Windows executable.

## Release contents

- Unified V1.0 application, executable, output-folder and build metadata.
- Full-covariance normal-reference distance, separate diagonal sensitivity distance, covariance entropy and distance-derived state scores.
- Advisory outlier detection without subject-level attenuation or exclusion; documented preprocessing and imputation still apply.
- Forest-module report wording aligned with modified Hartung-Knapp uncertainty intervals and endpoint classification.
- Bilingual supplementary methods distinguish the main workflow from the separate forest module, including their limitations and reproducibility requirements.
- Public source and archives contain templates and synthetic fixtures; research data and study results are excluded.

## Validation

The local Windows/Python 3.12 suite passed **138 tests**, with **7 optional external-data tests skipped**. Two existing Matplotlib layout warnings occurred in figure-export tests; those tests passed. Static checks passed. The source and packaged Windows executable both passed the synthetic end-to-end self-test (32 samples, 4 groups, 8 figures, 27 figure-data files), with exactly equal reference-distance summaries. Core numerical logic and the public numeric fixture were verified against the validated calculation baseline; changes to the forest module are limited to reporting and figure descriptions. The Windows file and product versions are 1.0.0.0.

The supplementary methods describe software implementation. Manuscripts must additionally specify the study's actual data, marker panels, reference groups, settings and comparisons. Software validation does not establish cross-model biological comparability.
