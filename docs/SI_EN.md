# BioEntropy V1.0 Supporting Information

Software methods and statistical definitions

## S1 Scope and input data

BioEntropy V1.0 provides normal-reference multivariate distance, covariance entropy, distance-derived state scores and secondary burden descriptors. A separate forest-plot module calculates standardized treatment contrasts within animal models and combines the resulting effect sizes. The two modules have different preprocessing and missing-data rules, which are specified separately below. These descriptors quantify features of the measured panel and do not by themselves establish recovery, causal treatment effects or physiological entropy.

In the main workflow, each workbook contains a group column, a sample identifier column and numeric biomarker columns. Longitudinal visits are supplied separately and identifiers can recur across visits. An animal workbook without a visit suffix receives the internal visit label 0W; this label does not establish a pretreatment sampling occasion. Treatment, comparator and normal-reference roles are selected for each run. A reference can comprise individual normal samples or uploaded lower and upper reference limits; interval limits do not supply a covariance matrix.

At each visit, an eligible main-workflow marker has at least 70% numeric nonmissing values and more than one distinct observed value. The intersection of eligible markers across visits defines the analysis panel, with at least two markers required. Eligibility is assessed across the records at a visit, not separately in each arm. The distance branch additionally requires a matched reference marker with at least three observed normal values. Let pE denote the anchor-preprocessing and covariance-entropy dimension, pN the matched distance dimension, n the group size at a visit, and m the normal-reference size. These dimensions need not be equal.

## S2 Anchor preprocessing

For covariance entropy and baseline-relative burden, parameters are fitted once using all groups at week 0 or the earliest available visit. A marker uses ln(1+x) when its observed anchor values are strictly positive and their absolute sample skewness exceeds 1; otherwise no logarithmic transformation is applied. Transformed values are clipped to the anchor q and 1−q quantiles, with default q=0.01. The same fitted limits apply to all visits. The scaled coordinate is

$$
z_{ij}=0.6745\frac{\widetilde t_{ij}-b_j}{a_j}\tag{S1}
$$

Here t̃ is the clipped transformed value, b is the median of the untrimmed transformed anchor values, and a is their median absolute deviation (MAD). If a<10⁻⁹, the implementation substitutes IQR/1.349, then the sample standard deviation, then 1 when necessary; the multiplier 0.6745 remains. For covariance entropy, remaining missing coordinates are filled with the group-and-visit column median. This preprocessing is distinct from normal-reference fitting in S3.

## S3 Normal reference distance

Normal-reference Mahalanobis deviation (NMD) is computed from a separate reference-fitting branch, rather than from Eq. S1. Marker names are matched using normalized names and recognized aliases. This branch selects ln(1+x) when the pooled observed study and matched-reference values are strictly positive and have absolute skewness greater than 1. Study and normal-reference values then share that transformation. Name matching and the glucose, lipid and HbA1c unit-conversion heuristics require verification against assay metadata; matching alone does not establish biological or measurement equivalence.

Within each transformed normal-reference marker, reference observations are clipped at the reference median ±4 MAD, with scale fallbacks as in S2 and an IQR fallback for invalid bounds. This clipping is applied to the reference used for fitting, not to study observations. Missing reference coordinates are filled with the clipped reference-column median. The reference mean μN and Ledoit–Wolf covariance ΣN are fitted to the completed transformed matrix. The covariance estimator shrinks the empirical covariance S toward a scaled identity [1]:

$$
\widehat\Sigma=(1-\lambda)S+\lambda\frac{\operatorname{tr}(S)}pI\tag{S2}
$$

The shrinkage intensity λ is estimated from the data. Covariance regularization stabilizes the matrix but is not an additional adjustment of individual samples and is not itself an outlier-robust estimator. In this branch, there is no further unit-variance standardization before covariance fitting. Units, the transformation and the retained panel therefore form part of the score definition. For study vector yi, the primary distance is

$$
D_i=\sqrt{(y_i-\widehat\mu_N)^T\widehat\Sigma_N^{-1}(y_i-\widehat\mu_N)}\tag{S3}
$$

Missing or nonfinite study coordinates are replaced with their fitted reference means. A record missing every retained coordinate can therefore receive distance zero; such a value is an imputation consequence, not observed normality. The primary distance has no ±8 coordinate cap. NMD and NMD_covariance_raw contain the same individual distance; group NMD is its arithmetic mean. Departures above or below the normal centre can increase D, so proximity to the reference is not a directional benefit for every marker.

The diagonal sensitivity score uses the same transformed values and reference means, with the observed clipped-reference sample SD sj; an invalid or near-zero SD is replaced by 1. It is

$$
D_{i,diag}=\sqrt{\sum_{j=1}^{p_N}\operatorname{clip}\!\left(\frac{y_{ij}-\widehat\mu_{N,j}}{s_j},-8,8\right)^2}\tag{S4}
$$

This comparison changes both covariance handling and capping; it does not isolate the effect of correlations. Exported χ² tail areas, and the scaled F approximation when m>pN+1, are descriptive reference-distribution approximations. Reference fitting, regularization and imputation prevent treating them as exactly calibrated treatment p values. The fitted-reference distance percentiles likewise do not establish a clinical normal range.

## S4 Entropy and distance derived state scores

For each group and visit, Gaussian covariance entropy is computed from the anchor-preprocessed, median-completed matrix using a group-specific Ledoit–Wolf covariance:

$$
H_{cov,g}=\frac12\ln|\widehat\Sigma_g|+\frac{p_E}2\ln(2\pi e)\tag{S5}
$$

Natural logarithms give nats in the specified transformed coordinates. This is the entropy of the fitted Gaussian covariance model, not a distribution-free estimate of biological disorder. Fewer than two records give an undefined estimate. If the covariance determinant is nonpositive numerically, the implementation sums logarithms of eigenvalues above 10⁻¹⁰; this is a numerical fallback. The fields Entropy_covariance, Entropy_raw and Entropy_robust_covariance carry the same covariance-entropy estimate. The field names do not imply that preprocessing is absent or that different estimators were fitted.

The software also computes a group state score from the second moment of distance about zero:

$$
H_{state,g}=\frac12\ln\!\left(\frac1{n_g}\sum_{i=1}^{n_g}D_i^2\right)+\frac12\ln(2\pi e)\tag{S6}
$$

This score requires at least two nonmissing distances and a positive finite mean square. It is not the differential entropy of the observed distance distribution. Because the mean square equals the squared mean plus the variance with divisor n, the score combines reference displacement and within-group spread. It is mathematically coupled to distance. Covariance entropy does not substitute distance into its formula; the two quantities must be named separately.

When Eq. S6 is defined, the generic field Entropy and its default intervals use Hstate and Entropy_method records normal_reference_state_entropy. Otherwise, Entropy falls back to Hcov and the method field records ledoitwolf_covariance_entropy. The explicit covariance and reference-state columns remain separate. Analyses must check this method field before comparing or correlating generic Entropy outputs. The generic interval does not describe Hcov when Hstate is selected.

## Outlier quality control

Quality control retains robust z, IQR, Isolation Forest, local outlier factor and minimum-covariance-determinant distance checks. The robust-z threshold is |z|>4; IQR fences use 1.5 IQR. Isolation Forest and local outlier factor require n≥5 and use contamination=min(0.05,max(0.01,5/n)); the latter uses min(20,max(2,floor(n/3))) neighbours. The covariance-determinant check requires n≥max(5,pE+1) and uses a χ²0.999,pE squared-distance threshold. An unavailable or failed check returns no flags. At least two flags define consensus and at least four define a strong flag. Flags are reported as quality-control counts and percentages; they do not delete samples or attenuate coordinates toward a group centre. Anchor clipping and missing-value handling still apply.

Diagnostics include n, pE, n/pE and shrinkage intensity. Labels reliable (n/pE≥5), limited (2≤n/pE<5) and low (n/pE<2) are implementation heuristics, not validated accuracy guarantees. Small n relative to p can make covariance entropy sensitive to regularization even when the matrix is invertible.

## S5 Secondary burden descriptors

The baseline-relative burden index (BRI, also exported as HDI) averages direction-adjusted coordinates from Eq. S1 over the pE markers. Missing coordinates contribute zero. The direction multiplier is −1 for recognized beneficial markers such as HDL-C and +1 for recognized adverse markers or unrecognized names. Non-HDL markers are checked before HDL to avoid substring misclassification. BRI is relative to the anchor population, not an absolute health scale. Direction rules require confirmation for the study context.

$$
BRI_i=\frac1{p_E}\sum_{j=1}^{p_E}s_jz_{ij},\quad s_j\in\{-1,+1\}\tag{S7}
$$

For supplied reference intervals, the normal-range burden score (NRBS) is the mean ln(1+relative excess) over observed supported markers; NRPS=1/(1+NRBS). For a two-sided interval [Lj,Uj], excess is zero inside the interval and the distance to the nearest bound outside it divided by max(Uj−Lj,0.1|Uj|,0.1|Lj|,10⁻⁹). With only an upper or lower bound, excess beyond that bound is divided by max(|bound|,10⁻⁹). Zero available markers produce undefined scores. Built-in adult intervals may be used when the panel matches the built-in fallback rules and no usable normal reference or uploaded interval exists; they require population and unit validation and are not animal normal ranges.

## S6 Uncertainty and graphical summaries

Nominal 95% bootstrap intervals use the 2.5th and 97.5th percentiles. Group means and selected Entropy use 200 resamples, treatment-minus-comparator mean differences use 400, entropy differences use 300, and the fitted normal-reference mean distance uses 500. The base random seed is 42 with deterministic derived seeds. Group resampling retains the original number of rows; contrast resampling samples the two groups separately. Reference means, covariance, transformations, clipping and initial imputations are held fixed. For Hcov, the group covariance is re-estimated within each resample of the preprocessed matrix. These intervals are conditional on the fitted preprocessing and reference, and do not propagate their full estimation uncertainty.

Burden contrasts also receive two-sided Mann–Whitney tests. For entropy differences, the exported bootstrap sign summary is twice the smaller proportion of resampled differences ≤0 or ≥0; this is a descriptive sign-proportion statistic rather than a separately calibrated null-resampling test. Spearman associations, visit-level linear trends and within-experiment standardized pooled associations are exploratory summaries. They are not subject-level repeated-measures models. Selected p-value families receive Benjamini–Hochberg adjustment [2], with adjusted columns identified by the fdr_bh suffix; this is not a single correction over every output.

Raw group distance plots display mean D; state-score plots display Eq. S6. Individual state-score dots can use ln(D)+½ln(2πe) for positive D, whose arithmetic mean generally differs from Eq. S6. Zero distances are omitted from these logarithmic overlays. Optional normalized displays anchor normal and model means at 0 and 100 when available, otherwise use within-experiment min–max scaling; equal values map to 50. The linear index L is clipped to [0,100], and contrast enhancement applies 100(L/100)³. These display transforms and separately scaled individual overlays are not effect sizes; raw estimates and intervals govern quantitative interpretation.

## S7 Independent cross model effect size analysis

The composite forest-plot module reads one row per animal and calculates within-model treatment-versus-PM contrasts before pooling across models. It is a separate calculation path: it does not import the main workflow NMD or Entropy values. In this module, Control denotes the normal reference and Model denotes the disease control. The implementation recognizes DB, KK-Ay and hamster datasets, with the medium KK-Ay dose selected for its default primary pool; separate low- and high-dose pools reuse the other model estimates. These overlapping dose pools are not independent replications.

Among recognized primary groups, markers with more than 20% missing values are removed. Complete-case rows are then selected across the remaining markers, before the subsequent main-panel direction filter. Thus rows are not recovered simply because a marker is later excluded. At least three complete normal-control animals are required. Markers with zero normal-control sample SD are dropped. Numeric entries bearing an asterisk are parsed and retained where possible, with their annotations recorded.

For the default composite, objective endpoints are retained only when the observed Model-versus-Control mean change matches the configured adverse direction. Context-dependent or unknown endpoints are excluded from the default composite. This is an observed-data inclusion rule, not a formal test that a model is valid, and it does not establish prospective endpoint prespecification. Signs are aligned consistently across groups. Each retained raw marker is centered and divided by the normal-control sample SD, and all resulting coordinates, including normal controls, are clipped to [−3,3]. This module does not apply the main workflow log-selection or anchor-quantile preprocessing.

Let u denote this clipped standardized vector. A Ledoit–Wolf covariance is estimated from the normal-control u vectors. The distance is computed as follows, using a matrix pseudoinverse if inversion fails:

$$
D_{i,F}=\sqrt{u_i^T\widehat\Sigma_{C,F}^{-1}u_i}\tag{S8}
$$

The distance centre is the zero vector defined by the original normal mean and SD. Although covariance fitting centers its input, the distance calculation does not subtract a newly estimated mean after clipping. A separate unrestricted sensitivity output removes the ±3 cap and uses all eligible endpoints without the main direction-inclusion filter. Because both panel selection and clipping differ, that output does not isolate either change alone.

The forest output labelled entropy is based on the per-animal log-distance score

$$
E_i=\ln\!\left(\max(D_{i,F},10^{-9})\right)+\frac12\ln(2\pi e)\tag{S9}
$$

Hedges’ g is calculated from E, not from group covariance entropy. The additive constant cancels in g. The module also exports descriptive group covariance entropies, but these do not enter this pooled entropy effect. The average of the individual scores in Eq. S9 generally differs from the group state score in Eq. S6. Linear-distance and log-distance effects summarize related information and cannot be treated as two independent confirmations.

## Effect size and pooling formulas

For either per-animal score x, a positive effect favors a smaller treatment score relative to PM. With treatment T, comparator C=PM, ν=nT+nC−2 and a nonzero pooled SD, the implementation uses the approximate Hedges correction [3]:

$$
s_p^2=\frac{(n_T-1)s_T^2+(n_C-1)s_C^2}{\nu},\quad d=\frac{\overline x_C-\overline x_T}{s_p}\tag{S10}
$$

$$
J=1-\frac3{4\nu-1},\quad g=Jd,\quad v=J^2\left(\frac{n_T+n_C}{n_Tn_C}+\frac{d^2}{2(n_T+n_C)}\right)\tag{S11}
$$

Nonfinite observations are excluded from each score comparison; groups require at least two finite values. The per-model interval is g±1.96√v. This variance approximation treats the derived animal scores as observed values; it does not explicitly propagate uncertainty from the fitted normal reference, endpoint selection or clipping.

Across k available models, fixed-effect inverse-variance weights are aj=1/vj. With fixed estimate gF, Q=Σaj(gj−gF)² and C=Σaj−Σaj²/Σaj, the DerSimonian–Laird estimate is τ²=max{0,[Q−(k−1)]/C} [4]. Random-effects weights and the pooled estimate are

$$
w_j=\frac1{v_j+\tau^2},\quad \widehat g_R=\frac{\sum_{j=1}^k w_jg_j}{\sum_{j=1}^k w_j}\tag{S12}
$$

A normal-approximation interval uses SE=1/√Σwj. For k≥2, the additional modified Hartung–Knapp interval uses q=Σwj(gj−ĝR)²/(k−1), SEHK=√[max(1,q)/Σwj], and ĝR±t0.975,k−1 SEHK [5]. For three models the implemented t multiplier is 4.303. The max(1,q) safeguard prevents a smaller adjusted standard error than the conventional random-effects value. I²=max{0,[Q−(k−1)]/Q}×100% when Q>0. With one contributing model, there is no cross-model heterogeneity information and no HK interval.

Different animal counts are accommodated through g, its variance and inverse-variance weighting; equal sample sizes are not required. Different marker panels are a separate limitation: standardization does not establish that distinct composites measure the same biological outcome. Cross-model pooling is exploratory unless that comparability is justified. Per-model panels, actual complete-case counts, dose choices and both pooled intervals should accompany the results. Common-marker and clipping sensitivity analyses can assess this dependence but are not automatically performed as a common-panel primary analysis. Correlated endpoints within a model, reused controls and overlapping dose pools must not be treated as independent studies.

## S8 Parameters and reproducibility

Table S1 summarizes implementation defaults rather than documenting the settings of a particular study. An article applying the software should report the actual retained panels, input units, normal-reference provenance and eligibility, sample counts after missing-data handling, selected groups and doses, visits, and any departures from these settings. Analysis records should retain input checksums, the source or executable checksum, release identifier, dependency versions and fitted parameter tables. The reference data and unit metadata determine the meaning of normal proximity.

**Table S1 Default calculation settings**

| Setting | Implementation default |
| --- | --- |
| Main marker coverage | ≥70% numeric nonmissing; >1 distinct value |
| Anchor quantile clipping | q=0.01; limits fitted once |
| Normal-reference fitting | Reference median ±4 MAD; scale fallbacks |
| Main distance | Full covariance; study coordinates uncapped |
| Diagonal sensitivity | Reference-SD coordinates capped at ±8 |
| Outlier voting | ≥2 methods consensus; ≥4 strong; flags only |
| Resampling and base seed | 200 / 400 / 300 / 500 as in S6; seed 42 |
| Forest missing-data rule | Drop markers >20% missing; then complete cases |
| Forest default distance | Normal-SD standardization; coordinate cap ±3 |
| Forest pooling | Hedges g; DL random effects; modified HK interval |

Lower D means proximity to the fitted reference; lower Hstate means a smaller distance second moment; lower Hcov means a smaller fitted covariance volume in the chosen coordinates. None alone defines biological recovery. Scores are not directly interchangeable across dimensions, panels, transformations or reference populations. Ledoit–Wolf regularization and bootstrap intervals do not remove small-sample bias, missing-data assumptions, dependence or multiple-comparison limitations. Source code and release materials are available at https://github.com/MrElss/BioEntropy. The exact release and checksum used for an analysis identify its implementation.

## References

1. Ledoit O, Wolf M. A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis. 2004;88:365–411. https://doi.org/10.1016/S0047-259X(03)00096-4

2. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

3. Hedges LV. Distribution theory for Glass’s estimator of effect size and related estimators. Journal of Educational Statistics. 1981;6:107–128. https://doi.org/10.3102/10769986006002107

4. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7:177–188. https://doi.org/10.1016/0197-2456(86)90046-2

5. Röver C, Knapp G, Friede T. Hartung–Knapp–Sidik–Jonkman approach and its modification for random-effects meta-analysis with few studies. BMC Medical Research Methodology. 2015;15:99. https://doi.org/10.1186/s12874-015-0091-1
