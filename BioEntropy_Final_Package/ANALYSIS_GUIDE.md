# BioEntropy V1.0 Analysis Guide

Before running an analysis, verify the input units, species, subject identifiers,
group roles and Normal-reference source. Keep a copy of the input workbooks and
record the selected parameters. The main workflow and the separate forest module
have different preprocessing and outcome definitions.

Multi-visit experiments output eight main artifacts. Figures 02 and 04
show overlapping nonzero starting points; duplicate figures 09/10 are omitted.
Single-visit experiments retain their nine artifacts and original figures 02/04.
For multi-visit displays, groups subtract their own baseline and add the mean of eligible group baselines, then
share a min/max 0–100 mapping and gamma=3 display. Decreasing or unchanged
trajectories start at 100. With a nonconstant range and no decrease, the lower
anchor is extended downward by the observed span; a constant range uses linear
50. Missing 0W uses the labelled earliest valid visit; fewer than two valid
visits do not define a trajectory. Burden dots share group offsets and anchors;
Entropy has no dots. The index is not a raw value or percentage, and zero does
not mean zero distance or entropy. It does not add paired or baseline-adjusted
inference. Raw figures 01/03 and result workbooks retain original estimates.

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
