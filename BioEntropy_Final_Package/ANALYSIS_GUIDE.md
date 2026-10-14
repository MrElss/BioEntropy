# BioEntropy V1.0 Analysis Guide

Before running an analysis, verify the input units, species, subject identifiers,
group roles and Normal-reference source. Keep a copy of the input workbooks and
record the selected parameters. The main workflow and the separate forest module
have different preprocessing and outcome definitions.

For longitudinal changes, use main figures 09 and 10. Their vertical axis is the
group metric minus its own 0W value, so eligible groups start together at zero.
If 0W is absent, use the recorded earliest valid visit; fewer than two valid visits
do not define a trajectory. These are changes in group summaries, not paired
individual changes or a statistical adjustment for baseline differences.
The raw plots retain observed between-group baseline differences.

Normalized subject points use the same group-summary anchors as the corresponding
curves or bars. Clipping and gamma=3 are display operations, and the transformed
group value need not equal the mean of transformed points. Entropy plots show
only group estimates, without individual log-distance dots.

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
