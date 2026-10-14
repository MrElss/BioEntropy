# BioEntropy V1.0 — display correction

The application remains **BioEntropy V1.0**. Release **v1.0.1** and Windows file
version **1.0.1.0** identify this corrected build.

## What changed

- Longitudinal outputs now include **Figure 09: disease-deviation change** and
  **Figure 10: Entropy change**. Each group subtracts its own 0W estimate, so the
  baseline markers coincide at zero. If 0W is unavailable, the earliest valid
  visit is used and recorded. Groups with fewer than two valid visits do not
  define a change trajectory. Single-timepoint outputs remain group comparisons.
- Individual disease-deviation points now use the same group-summary reference
  anchors as the corresponding curves or bars. The original normalized group
  values are preserved. Clipping and gamma=3 mean a transformed group value
  still need not equal the mean of the transformed individual points.
- Entropy figures show group estimates and applicable existing intervals;
  individual log-distance dots are no longer overlaid as entropy observations.
- Figure-data CSV files follow the actual displayed values and include the
  baseline week, baseline value and change. Files for unplotted point layers
  are no longer generated. New baseline-figure labels remain readable at crossings.
- Chinese and English supplementary methods, input guidance and release metadata
  describe the corrected display rules.

## What did not change

The NMD, covariance entropy, reference-state score, NRBS/BRI, outlier/QC,
bootstrap and forest-module calculations are unchanged. Raw numerical tables,
existing statistical results and original normalized group curves are preserved.
The new baseline figures describe changes in group summaries; they do not add a
paired-subject analysis, a new uncertainty interval or a baseline-adjusted
treatment-effect estimate. The original absolute-value plots remain available.

## Validation

- **155 tests passed; 7 optional external-data tests skipped.** Static checks passed.
  Two existing Matplotlib layout warnings occurred in unrelated figure-export tests.
- Numerical fixtures are unchanged, and the public regression suite checks the
  full result tables against those fixtures. Computational modules were also
  compared with the preceding release.
- Source and Windows executable end-to-end self-tests passed: 32 synthetic
  sample records, 4 group/visit summaries, 10 figures and 21 figure-data files.
  Raw distance and entropy summaries agree between source and executable.
- Existing longitudinal study tables were checked locally: baseline changes
  begin at zero, normalized points use the group-curve anchors, original group
  values match and the study files remain untouched. Those data are not distributed.
- Bilingual Word documents were checked for method consistency and page layout.

## Downloads

- **BioEntropy_V1.0_Windows.zip** — extract and run `BioEntropy_V1_0.exe` on Windows x64.
- **BioEntropy_V1.0_Source.zip** — source, tests, templates and documentation.
- **BioEntropy_V1.0_SI_CN.docx / BioEntropy_V1.0_SI_EN.docx** — editable supplementary methods.
- **BioEntropy_V1.0_SHA256.txt** — checksums for the downloads and executable.

For overlapping 0W starting points, open **main figure 09 or 10**. The original
absolute-value and normalized figures retain observed baseline differences.
