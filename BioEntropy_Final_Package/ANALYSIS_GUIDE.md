# BioEntropy V1.0 Analysis Guide

Before running an analysis, verify the input units, species, subject identifiers,
group roles and Normal-reference source. Keep a copy of the input workbooks and
record the selected parameters. The main workflow and the separate forest module
have different preprocessing and outcome definitions.

The primary NMD uses a Ledoit-Wolf full-covariance reference. The diagonal score
with ±8 standardized-value capping is a sensitivity output. Outlier flags are
quality-control information and do not trigger subject-level attenuation or exclusion.
Other preprocessing, including winsorization and missing-value imputation, remains active.

Always state whether an entropy result is group covariance entropy, the group
distance-derived state score, or the forest module's per-animal log-distance score.
Consult `Entropy_method` and the metric dictionary. Distance-derived scores are
not independent of NMD. Cross-model composite pooling is exploratory and requires
an explanation of marker-panel comparability.

Report study-specific sample sizes, included markers, missing-data handling,
Normal-reference selection, group contrasts, bootstrap settings and sensitivity
analyses. Software defaults alone do not establish which settings a study used.
The repository `docs/` directory contains Chinese and English supplementary methods.
