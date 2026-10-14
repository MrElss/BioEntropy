# BioEntropy V1.0 — Methods and Formulas

This technical reference describes the indices implemented in BioEntropy and
their statistical treatment. The source and run audit records identify the
implementation and parameters used for an analysis. Notation: a group `g` has `n` samples on `p` common numeric
features; `Σ_g` is the group covariance; `μ̂, Σ̂` are estimated from a Normal
reference of `m` samples.

## 1. System entropy (log-determinant Gaussian entropy)

For `X ~ N(μ, Σ)` in `p` dimensions the differential entropy is

```
H(X) = ½ · log[ (2πe)^p · |Σ| ] = ½·log|Σ| + (p/2)·log(2πe)
```

We estimate `Σ` per group with a **Ledoit–Wolf shrinkage** covariance (so it is
invertible even when `n` is small) and report `H_g`. Within one experiment the
constant `(p/2)·log(2πe)` cancels in group comparisons.

**Reliability (n/p).** The estimator needs `Σ̂_g` to be well conditioned. When
`n ≤ p` the sample covariance is singular and `log|Σ̂_g|` is dominated by the
shrinkage *target*, not the data. We therefore report, per group:
`Entropy_cov_n_samples (n)`, `Entropy_cov_feature_count (p)`,
`Entropy_cov_n_over_p`, the Ledoit–Wolf shrinkage intensity, and a label:

| n/p        | label      | interpretation                                   |
|------------|------------|--------------------------------------------------|
| ≥ 5        | `reliable` | covariance well estimated                        |
| 2 – 5      | `limited`  | usable but shrinkage-influenced                  |
| < 2        | `low`      | covariance singular/over-shrunk; compare with care |

When a true Normal reference exists, the pipeline reports the
**normal-reference state entropy**, an entropy-like score using the second
moment of NMD around zero:

```
H_state = ½·log(mean(NMD²)) + ½·log(2πe)
```

It is not Gaussian entropy calculated from the variance of the NMD distribution:
shifting all distances while preserving variance changes this score. It avoids
estimating a p-dimensional covariance per group, but still depends on the fitted
Normal reference and missing-value treatment. The primary distance is not z-capped.

> **Circularity caveat.** The reference-state score directly reuses each group's
> mean squared NMD, and is monotonic in that second moment. Associations with
> NMD therefore have structural coupling; it is not the differential entropy of
> the NMD distribution. A positive
> entropy–deviation association computed with reference-state entropy is
> therefore **not independent corroboration** of the hypothesis. The
> **covariance** entropy (`Entropy_robust_covariance`) does not directly reuse
> NMD in its formula; use it for a separately defined association, and report the
> n/p reliability alongside it. `tests/test_methodology_validation.py` pins both
> behaviours on synthetic data.

### 1a. Subject handling in V1.0

Baseline transformation, default 1%/99% winsorisation and median imputation remain.
No individual is moved toward a median or excluded because of outlier votes.
The five methods still record >=2 consensus and >=4 strong QC flags. Flags do
not affect covariance fitting, entropy bootstrap inputs or sample inclusion.
Ledoit-Wolf regularizes a covariance matrix; it is not subject attenuation and
does not automatically provide robustness to outliers. `Entropy_raw`,
`Entropy_covariance` and the legacy alias `Entropy_robust_covariance` are equal.
`Entropy_n_adjusted_subjects` and legacy attenuation-count fields are zero;
`Entropy_n_strong_qc_flags` reports flagged subjects separately.

## 2. NMD — Normal-referenced Mahalanobis deviation

Primary NMD uses the fitted **Ledoit-Wolf full Normal-reference covariance**:

```
NMD² = (x − μ̂)ᵀ Σ̂_LW⁻¹ (x − μ̂)
```

Normal-reference matching, unit harmonization, optional log1p, Normal-only
median ± 4 MAD clipping and missing-value imputation are applied to the Normal reference.
Normal clipping fits the reference; study observations are not clipped to these
bounds. The primary distance uses transformed study observations with no ±8 z cap.
`NMD_covariance_raw` is a compatibility alias of primary NMD.

The old diagonal score is retained as **`NMD_diagonal_capped`**: deviations from
the fitted Normal mean are divided by the sample SD of the processed Normal
observations, capped at ±8 and combined in quadrature. `NMD方法敏感性` reports
group means and bootstrap intervals; `治疗比较` includes diagonal mean differences,
intervals and exploratory Mann-Whitney p-values. These contrasts compare both
covariance handling and capping, so they do not isolate a covariance-only effect.

**Reference tail areas.** `NMD_p_chi2 = chi2.sf(NMD², p)` now describes the full
distance. `NMD_diagonal_p_chi2` retains the previous diagonal reference output.
`NMD_p_smallsample_F` aliases `NMD_cov_p_smallsample_F`, using the existing
predictive F scaling only for m > p + 1. Neither chi-square nor F calibration is
exact with estimated, clipped/imputed, Ledoit-Wolf reference covariance. These
outputs are exploratory, and are not the treatment-comparison p-values.

### 2b. Normal-reference NMD band (reported, not shaded)

`compute_normal_reference_nmd_band` reports the distance-from-Normal that the
healthy reference subjects themselves carry (the same full-covariance Mahalanobis
transform applied to the Normal cloud): healthy per-subject central 95% range
(p2.5–p97.5) and the group mean with a bootstrap CI. It describes distances in the fitted Normal reference and does not establish
statistical equivalence to Normal. The fitted-reference summaries are
reported with the weekly state display rather than used as an equivalence threshold.

### 2c. Weekly state small-multiples

An alternative clinical trajectory view (main figure 08): one panel per week, with
each arm placed as a bubble in the (entropy, burden) state plane on axes shared
across all panels, so the week-to-week march is directly comparable. Diseased
corner is top-right, healthy corner bottom-left; a thin connector plus a
significance star (burden-gap Mann–Whitney p) marks the per-week separation. This
makes the "treatment pulls toward healthy first, then the crossed-over comparator
catches up" pattern legible at a glance.

### 2a. Ingesting an uploaded Normal-sample file

A Normal cloud can be supplied as a separate workbook (e.g. `Normal.xlsx`) rather
than an in-file group. The loader accepts:

* a **prefix-less** `Normal.xlsx` (any `TRUE_NORMAL_FILE_ALIASES` stem) as a
  *global* reference applied to every experiment in the batch, or an
  experiment-scoped `105-Normal.xlsx`;
* an **ID-only layout** (an ID column followed by marker columns, no Group
  column) — parsed as one `Normal` cloud;
* **unit mismatches** — glucose/lipid markers delivered in conventional US units
  (mg/dL) are auto-converted to SI (mmol/L) via known factors (glucose ×0.05551,
  cholesterol ×0.02586, triglyceride ×0.01129). The conversion is adopted only
  when same-unit medians are biologically implausible and the converted median
  lands within a plausible ratio of the study median; the factor and a note are
  recorded in the Normal-reference parameter table.

A Normal file whose markers do not overlap an experiment's panel is skipped for
that experiment (it falls back to §3a/§3), never aborting the run.

## 3. BRI / HDI — baseline-relative burden

Direction-corrected, robust (median/MAD) z-scores aggregated across features and
anchored to the earliest week. `infer_direction(feature)` sets the sign so that
"higher = worse" markers and "higher = better" markers (e.g. HDL-C) point the
same way. BRI is a **relative** burden index, not an absolute health score.

Because BRI is anchored to each cohort's own baseline week (no healthy floor), it
can misframe a placebo/untreated arm: shared regression-to-the-mean and
standard-of-care drift render as apparent "improvement", and small absolute
differences near the low end can cross over. Prefer an absolute anchor (§3a/§2)
when one is available.

### 3a. NRBS — Normal-range burden score (absolute anchor)

When a study has **no enrolled Normal/healthy arm** and **no uploaded
Normal-range file**, clinical-lab markers are anchored to built-in **published
adult reference ranges** (`PUBLISHED_LAB_REFERENCE_RANGES`), giving an *absolute*
health scale instead of a baseline-relative one:

```
excess_j = relative distance of marker j outside its published Normal range
           (0 when in range; one-sided — only the pathological bound is set)
NRBS     = mean_j log1p(excess_j)          # 0 ⇒ all markers within Normal range
```

Direction is encoded by which bound is open: higher-worse labs (FPG, HbA1c, TG,
LDL-C, non-HDL-C, TC, ALT, AST, GGT, hsCRP) set an **upper** limit; HDL-C sets a
**lower** limit. Marker→range matching uses the longest (most specific) token, so
`非高密度脂蛋白胆固醇` binds to the non-HDL limit, not the HDL limit. Panels with
no matching lab (e.g. animal studies) get an empty table and fall back to BRI.

**Metric preference for the main figures** (`_preferred_burden_kind`):
`NMD` (true enrolled Normal cloud) > `NRBS` (published-range anchor) > `BRI`
(baseline-relative fallback).

### 3b. Treatment-minus-comparator contrast (results tables)

Δ = (treatment − comparator) burden per week with a nominal 95% bootstrap CI.
Δ < 0 means a lower observed group-average burden in the treatment arm.
Causal interpretation depends on the study design and assumptions; this
contrast alone does not adjust for baseline differences or identify a drug effect.
The results also report per-week Mann–Whitney comparisons.

## 4. Effect size — Hedges' g (forest module)

```
pooled SD = sqrt(((n1−1)·SD1² + (n2−1)·SD2²) / (n1+n2−2))
Cohen's d = (mean1 − mean2) / pooled SD
J         = 1 − 3 / (4·(n1+n2−2) − 1)
Hedges' g = J · Cohen's d
```

## 5. Association testing and multiplicity

The core hypothesis (entropy ↔ disease deviation) is summarised by Spearman
correlations. These are computed **across groups**, so the effective sample size
is the number of groups (often 4–7) — treat them as **exploratory/descriptive**,
not confirmatory.

**FDR.** p-values are adjusted with **Benjamini–Hochberg** within two
pre-specified families and reported as `<col>_fdr_bh` q-values:

1. deviation-vs-entropy associations (`实验内相关性` + `跨实验合并相关性`);
2. metric-vs-time trends (`时间趋势`).

HDI mirrors BRI and is excluded from the family to avoid double-counting.

## 6. Reproducibility

Bootstrap and stochastic quality-control methods use a fixed seed (`RNG_SEED = 42`).
`tests/test_public_pipeline.py` checks fixed synthetic inputs against
`tests/public_numeric_baseline.json`, including numeric tables, workbook exports,
figure generation, output layout and audit isolation. `test_calculation_methods.py`
checks covariance orientation, distance formulas, extreme values, missing inputs,
regularization and the absence of subject-level attenuation. The independent
`reference_metric_baseline.json` fixture checks the primary and sensitivity distance
routes and the BRI, NRBS and covariance-entropy calculations.

Reproducibility requires the same inputs, settings, software and dependency versions.
Rendering also depends on the installed fonts. Public synthetic checks do not replace
study-specific validation. Run manifests record the version, method policy, source
and input hashes, and available environment information.

## 7. Interpretation

### Display and baseline changes

Longitudinal figures 02/04 align group baselines: A(g,t)=X(g,t)−X(g,b)+B,
where B is the arithmetic mean of eligible group baselines. Prefer a finite 0W
estimate; otherwise label the earliest valid visit. Groups require two valid
visits. Aligned group values fit common min/max anchors within each experiment
and metric. Clip the linear 0–100 index L and display 100*(L/100)^3.
If no change is negative and the range is nonconstant, extend the lower anchor
downward by the observed span so the baseline stays nonzero. A constant range
uses L=50. With decreasing or unchanged trajectories, the shared baseline is 100.
Burden dots share their group's offset and the curve's anchors; Entropy has no
individual dots. Raw estimates, differences and anchors are exported. This is
not a raw metric or percent change; zero does not mean zero distance or entropy.
No paired inference, new intervals or baseline-adjusted treatment effects are added.
Raw metrics and existing statistical outputs remain unchanged.

For normalized level figures, the group curves/bars and disease-deviation subject
points use the same anchors fitted to that figure's valid group summaries within
the experiment: distinct Normal/model means when available, otherwise their
minimum/maximum. Joint views can first require both metrics to be present and
therefore use a smaller subset than a single-metric view.
The common linear mapping is clipped to [0,100] and optionally transformed by
100*(L/100)^3. Single-timepoint normalized state figure 08 uses linear L without
gamma. A constant reference range gives L=50 (12.5 after the cubic transform).
Thus a transformed group mean is not generally the mean of the
transformed points; endpoint values may reflect clipping rather than observed
normality. Entropy remains a group-level estimate and is displayed without
individual log-distance points. The normalized group curves retain their prior
values; the baseline-aligned normalized figures also use gamma=3.

### Scientific interpretation

The distance-derived state score is mathematically coupled to NMD. Their correlation
does not constitute independent validation. Forest-module log-distance effect sizes
are distinct from group covariance entropy. Cross-model pooling requires comparable
outcome definitions; standardization alone does not establish this comparability.
With few models, report the model-specific estimates and modified Hartung-Knapp
intervals and treat composite pooled results as exploratory. The standalone forest
module uses its own preprocessing and must be described separately from the main
analysis. See the bilingual supplementary methods in `docs/` at the repository root.

Multi-visit experiments export eight main artifacts; baseline-aligned figures 02/04 replace the former duplicate 09/10 outputs. Single-visit experiments retain nine main artifacts and original normalized level figures 02/04. Raw estimates remain in figures 01/03 and result workbooks.
