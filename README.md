# BioEntropy V1.0.2

BioEntropy is a Windows application for multivariate analysis of clinical and animal data. It provides Normal-reference Mahalanobis distance (NMD), Shannon entropy, Rényi entropy of order 2, energy distance, Hedges’ g and forest plots.

## Windows application

Extract **BioEntropy_V1.0.2_Windows.zip** and run **BioEntropy_V1.0.2.exe**.

Load Excel files or folders. The application detects grouped and longitudinal experiments automatically. Set the Normal reference, treatment, comparator and display order, then select a method and start analysis. The result folder defaults to the Desktop; repeated names receive a numeric suffix. Chinese and English are available for the interface and output text.

Each analysis provides a results workbook, figures, plotting workbooks, saved settings, a run log and an HTML results viewer. Plotting workbooks contain group estimates, sample scores, Prism tables and P values by figure. Source labels and patient identifiers retain their original values.

## Methods

- **NMD:** multivariate deviation from the Normal reference.
- **Shannon entropy:** uncertainty of the smoothed joint marker distribution at a fixed common scale; scatter points show sample information scores.
- **Rényi entropy of order 2:** a complementary density-based measure using the same fixed-scale kernel distribution.
- **Energy distance:** distributional discrepancy from Normal, incorporating location, spread and shape.
- **Hedges’ g:** standardized effects and forest plots.

Select Shannon, Rényi, energy distance or NMD only. The selected method determines the output. Statistical settings, inclusion counts and P values are recorded in the results workbook.

With a usable Normal reference, Shannon and energy distance use a common eligible cohort: complete Shannon-panel records with valid NMD, with at least two records per group and visit. Identical input and reference settings yield identical NMD sample values, counts and group means in these two modes. The results workbook contains a common-sample ledger with identifiers, inclusion decisions, reasons, missing markers and source rows. NMD Holm-adjusted P values may differ across methods because multiplicity adjustment includes the selected metric. Rényi uses complete entropy-panel records. NMD-only and the independent G-value module apply their respective eligibility rules.

## Source

Python 3.12 on Windows:

```text
python -m pip install -r requirements.txt
python BioEntropy.py
```

Command-line analysis:

```text
python BioEntropy.py --run --input "input_folder" --output "new_result_folder" --display Shannon --permutations 9999 --language en
```

`--display` accepts `Shannon`, `Renyi2`, `Energy` and `none`. `--project` loads saved analysis settings. `--language` accepts `zh` and `en`.

## Tests and build

```text
python -m pip install -r requirements-dev.txt
python -m ruff format --check .
python -m unittest discover -s tests -v
python -m PyInstaller --clean --noconfirm BioEntropy.spec
```

## 中文使用

下载并解压Windows软件包，运行BioEntropy_V1.0.2.exe。载入Excel文件或文件夹，确认正常参照、治疗组、比较组、显示顺序和计算方法后开始分析。软件自动识别实验类型，结果默认保存到桌面，重名目录自动追加序号。详细操作见[使用说明](使用说明.md)，方法和散点定义可在软件内查看。

## License

See [LICENSE](LICENSE), [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and the `licenses` directory.
