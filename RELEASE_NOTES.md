# BioEntropy V1.0 figure consolidation

The application remains **BioEntropy V1.0**. Release **v1.0.3** and Windows file version **1.0.3.0** identify this build.

## Output changes

- Multi-visit experiments now use the baseline-aligned normalized NMD/NRBS/BRI display as **main figure 02**, and the baseline-aligned normalized Entropy display as **main figure 04**. These are the former figures 09 and 10; their displayed values and scaling rules are unchanged from v1.0.2.
- Duplicate longitudinal figures 09/10 and their duplicate figure-data exports are removed. Multi-visit experiments ordinarily produce **eight main artifacts**. Raw figures 01/03 and raw workbooks retain the original estimates.
- Single-visit experiments retain **nine main artifacts**, including their original normalized level figures 02/04 and figure 09 for mechanism drivers. No longitudinal alignment is applied to them.
- Figure-data filenames, element mappings, desktop guidance and bilingual supplementary methods follow the new numbering. Individual burden dots share the group baseline offset and display anchors; Entropy is shown by group estimates only.

## Validation

Scientific formulas, raw result tables, intervals, comparisons, effect sizes and other figures remain unchanged. The aligned display index is neither a raw metric nor a percentage change.

- **162 tests passed; 7 optional external-data tests skipped.** Static checks passed. Two existing Matplotlib layout warnings remain in unrelated single-visit export tests.
- Source and packaged Windows self-tests passed: 32 synthetic sample records, 4 group/visit summaries, 8 main artifacts and 19 figure-data files, with exactly equal raw distance and entropy summaries.
- Experiment 106 was rerun locally against the supplied preceding result archive. All 10 workbooks, including 142,108 numerical cells and their text fields, match exactly. New figures 02/04 match the preceding figures 09/10 in every exported curve value; all 4,700 NMD display points also match. Inputs and the previous archive are untouched and are not distributed.
- Single-visit regression results and artifact layouts are preserved. Both supplementary Word documents have 8 pages and 13 editable equations; all equations are unchanged.

## Downloads

- **BioEntropy_V1.0_Windows.zip** — extract and run `BioEntropy_V1_0.exe` on Windows x64.
- **BioEntropy_V1.0_Source.zip** — source, tests, templates and documentation.
- **BioEntropy_V1.0_SI_CN.docx / BioEntropy_V1.0_SI_EN.docx** — supplementary methods.
- **BioEntropy_V1.0_SHA256.txt** — download and executable checksums.

For aligned nonzero longitudinal starting points, open **main figures 02 and 04**.
