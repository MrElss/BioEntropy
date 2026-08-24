# BioEntropy — Methods and Formulas (authoritative reference)

This file is the single source of truth for the indices BioEntropy computes and
their statistical treatment. Code comments and the auto-generated report should
defer to it. Notation: a group `g` has `n` samples on `p` common numeric
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
**normal-reference state entropy** (entropy of the 1-D NMD distribution), which
does not suffer the n/p problem.

> **Circularity caveat.** The reference-state entropy is the entropy of each
> group's NMD values, and groups further from Normal have NMD distributions that
> are both larger and more spread out. So the reference-state entropy correlates
> with disease deviation *almost by construction*. A positive
> entropy–deviation association computed with reference-state entropy is
> therefore **not independent corroboration** of the hypothesis. The
> **covariance** entropy (`Entropy_robust_covariance`) is the measure that is
> independent of deviation; use it for the headline association, and report the
> n/p reliability alongside it. `tests/test_methodology_validation.py` pins both
> behaviours on synthetic data.

## 2. NMD — Normal-referenced Mahalanobis deviation

Distance of a sample from the robust Normal envelope. Two variants are stored:

- **Diagonal (primary `NMD`)**: per-feature robust z-scores, capped at
  `±NMD_Z_CAP`, summed in quadrature; `NMD² = Σ_j z_j²`.
- **Full covariance (`NMD_covariance_raw`)**: `D² = (x-μ̂)ᵀ Σ̂⁻¹ (x-μ̂)`.

**p-values.**
- `NMD_p_chi2 = 1 − F_{χ²(p)}(NMD²)` — assumes the scales are *known* and
  features independent-normal; it is **optimistic** because μ̂/scale are
  estimated from a finite reference.
- `NMD_cov_p_smallsample_F` — exact small-sample p-value for the full-covariance
  distance. For a new observation versus a Normal reference of `m` samples,

  ```
  ((m − p) / (p (m − 1))) · (m / (m + 1)) · D²  ~  F(p, m − p)
  ```

  (out-of-sample Hotelling T²). It is defined only when `m > p + 1`; otherwise
  it is `NaN`, which itself signals that the covariance-based p-value is not
  trustworthy for that reference size.

### 2b. Normal-reference NMD band (reported, not shaded)

`compute_normal_reference_nmd_band` reports the distance-from-Normal that the
healthy reference subjects themselves carry (the same robust-capped Mahalanobis
transform applied to the Normal cloud): healthy per-subject central 95% range
(p2.5–p97.5) and the group mean with a bootstrap CI. It is the range a treatment
arm would have to reach to be indistinguishable from Normal. On chronic-disease
lab panels the treated arms typically stay well above it (e.g. exp 105 ≈ 7.7–9.6
vs a healthy band of ≈ 0.9–4.0), so it is surfaced as a compact note on the
weekly state small-multiples (§2c) rather than shaded on the main NMD trajectory,
where an all-empty band far below the data would only compress the axis.

### 2c. Weekly state small-multiples

An alternative clinical trajectory view (main figure 11): one panel per week, with
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

### 3b. Treatment-minus-comparator contrast (main figure 12)

Δ = (treatment − comparator) burden per week with 95% CI. The difference cancels
the shared placebo/standard-of-care drift, isolating the **drug-attributable**
effect (the controlled clinical estimand). Δ < 0 ⇒ the treatment arm carries
lower disease burden; per-week Mann–Whitney significance is starred.

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

All bootstrap/resampling uses a fixed seed (`RNG_SEED = 42`). The full pipeline
output (tables, reports, and every figure's pixels) is pinned by
`tests/test_regression_pipeline.py` against `tests/baseline_snapshot.json`;
regenerate the baseline only for an *intended* methodological change and review
the diff.

## 7. Caveats (publication)

Exploratory visual outputs are not a substitute for a predefined statistical
analysis plan. Define one primary endpoint (e.g. a subject-level entropy–NMD
association with experiment as a random effect), fix it in advance, and treat
all other outputs as exploratory.
