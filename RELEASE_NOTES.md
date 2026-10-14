# BioEntropy V1.0 — aligned nonzero baselines

The application remains **BioEntropy V1.0**. Release **v1.0.2** and Windows file version **1.0.2.0** identify this build. It supersedes the longitudinal display in v1.0.1.

## What changed

- Main figures 09 and 10 now align each group's baseline estimate to one common, nonzero point. For trajectories that decrease after baseline, both curves start at 100 on the display index. A valid 0W estimate is preferred; otherwise the earliest valid visit is identified in the legend. Two valid visits per group are required.
- Aligned group estimates use one shared range for each experiment and metric, are clipped to a linear 0–100 index, and then use the same gamma=3 display transform. Constant trajectories and trajectories without a decrease retain a nonzero shared start under the documented anchor rules.
- Individual NMD, NRBS or BRI dots use their group's baseline offset and the same range, clipping and transform as the group curves. Entropy is displayed by group estimates only. Exported figure data include raw estimates, raw baseline differences, aligned estimates and display anchors; the index is neither a raw distance/entropy value nor a percentage change.
- The Chinese and English supplementary methods and the software guidance describe the revised figures and their interpretation.

## Calculation and validation

Scientific calculations, existing group estimates, intervals, comparisons, effect sizes and original figures 01–08 are unchanged. The new figures are descriptive summaries and do not imply paired-subject inference or baseline-adjusted ANCOVA.

- **162 tests passed; 7 optional external-data tests skipped.** Two existing Matplotlib layout warnings remain in unrelated figure-export tests.
- Source and Windows executable self-tests passed using 32 synthetic sample records, 4 group/visit summaries, 10 figures and 23 figure-data files. Raw distance and entropy summaries agree exactly.
- Experiment 105 was rerun locally with the new source. All 10 result workbooks match the preceding run at every one of 100,693 numerical cells; maximum absolute difference is zero. Its new curves also match the accepted display preview. Private study inputs and outputs are not included in the download.
- Both supplementary Word documents contain 8 pages and 13 editable equations; the 12 original scientific equations are unchanged.

## Downloads

- **BioEntropy_V1.0_Windows.zip** — extract and run `BioEntropy_V1_0.exe` on Windows x64.
- **BioEntropy_V1.0_Source.zip** — source, tests, templates and documentation.
- **BioEntropy_V1.0_SI_CN.docx / BioEntropy_V1.0_SI_EN.docx** — editable supplementary methods.
- **BioEntropy_V1.0_SHA256.txt** — checksums for the downloads and executable.

For coincident nonzero baseline points, open main figure **09** or **10**. Figures 01–08 retain their original definitions.
