from constants import SOFTWARE_VERSION


METHOD_TITLES = {
    "zh": {
        "Shannon": "Shannon熵",
        "Renyi2": "Rényi二阶熵",
        "Energy": "能量距离",
        "none": "仅NMD",
    },
    "en": {
        "Shannon": "Shannon entropy",
        "Renyi2": "Rényi entropy of order 2",
        "Energy": "Energy distance",
        "none": "NMD only",
    },
}

COMMON_METHOD_GUIDE = {
    "zh": """NMD按正常参照的协方差计算多指标距离，组NMD为有效样本距离的平均值。

结果表.xlsx包含组结果、样本结果、原值与显示值、显著性、纳入人数、运行设置及输入来源。
图表目录保存结果图；绘图数据目录按实验和图号保存组值、散点、Prism数据和P值。
查看结果.html用于浏览结果。

图01/03显示原值，图02/04使用线性显示尺度，图05对照两种指标。
动物图07/08显示组状态；临床图10/11显示分周原值联合状态，保留0W原始坐标。
小图展示散点的完整分布。

临床图01/03和10/11比较当周原值，图02/04/05比较相对0W的变化；动物实验比较原值。
置换单位为整只动物或同一试验内的完整患者轨迹，采用Holm校正。
各图P值及适用条件记录在结果表中。数值积分误差反映计算精度。

结果表按组和访视列出各指标的纳入人数。
Shannon与能量距离使用相同的共同样本规则。共同样本清单逐条列明个体编号、纳入状态、原因、缺失指标及原表行号。
独立G值功能按全部可计算NMD的动物估计NMD效应，各组人数见其结果表。
""",
    "en": """NMD uses the Normal-reference covariance to measure multivariate deviation.
Group NMD is the mean of valid sample distances.

结果表.xlsx contains estimates, sample values, raw and display values, significance
tests, inclusion counts, settings and input sources. 图表 contains figures;
绘图数据 contains plotting workbooks by experiment and figure, including group values,
scatter scores, Prism data and P values. Open 查看结果.html to browse the results.

Figures 01/03 show raw values; figures 02/04 use a linear display scale;
figure 05 compares the metrics. Animal figures 07/08 show group states.
Clinical figures 10/11 show weekly joint states with the original 0W coordinates.
Insets show the full scatter distributions.

Clinical figures 01/03 and 10/11 compare raw values at each visit;
figures 02/04/05 compare changes from 0W. Animal tests compare raw values.
Permutation tests exchange whole animals or complete patient trajectories within
each trial and use Holm adjustment. P values and interpretation conditions are
saved in the results workbook. Numerical integration errors measure precision.

The results workbook lists sample counts by group, visit and metric.
Shannon and energy distance share the same sample eligibility rule. The common-sample ledger records each identifier, inclusion decision, reason, missing markers and source row.
The independent G-value analysis uses all animals with a valid NMD for NMD effects;
its workbook reports the corresponding group sizes.
""",
}

METHOD_DEFINITIONS = {
    "Shannon": {
        "zh": """H = −∫f̂(z)ln f̂(z)dz，单位为nats。
f̂为固定共同尺度和带宽的等权高斯核密度，描述多指标联合分布的不确定性。

有可用正常参照时，NMD与Shannon按组和访视使用同一批动物或患者记录，要求NMD有效、熵面板完整且组内至少两条有效记录。组均值、散点和置换检验采用该集合。
散点为样本信息分数：C_i = −∫K_h(z−z_i)ln f̂(z)dz。
C_i的平均值在数学上等于组Shannon熵，C_i/n为等权贡献；数值积分存在计算误差。
显著性检验在每次置换中重新计算组统计量。
""",
        "en": """H = −∫f̂(z)ln f̂(z)dz, in nats.
f̂ is an equally weighted Gaussian kernel density with a fixed common scale and
bandwidth. It describes uncertainty in the joint marker distribution.

With a usable Normal reference, NMD and Shannon use the same records per group and visit: a valid distance, a complete entropy panel and at least two eligible records per group and visit. Group estimates, scatter points and tests use this cohort.
Scatter scores: C_i = −∫K_h(z−z_i)ln f̂(z)dz.
The mean of C_i equals the group Shannon entropy mathematically; C_i/n is the
equally weighted contribution. Numerical integration introduces computational
error. Significance tests recompute the group statistic for each permutation.
""",
    },
    "Renyi2": {
        "zh": """H_2 = −ln∫f̂(z)²dz，单位为nats。
采用同一固定尺度核密度，通过高斯核重叠解析计算，侧重分布中的高密度区域。
图中展示组估计，临床图10展示分周组点。
""",
        "en": """H_2 = −ln∫f̂(z)²dz, in nats.
It uses the same fixed-scale kernel density, calculated analytically from Gaussian
kernel overlaps, and emphasizes high-density regions of the distribution.
Plots show group estimates; clinical figure 10 displays these estimates by visit.
""",
    },
    "Energy": {
        "zh": """E_U = 2mean||X−Y|| − Umean||X−X′|| − Umean||Y−Y′||。
Y为Normal参照，组内项使用不同样本之间的配对，原值允许为负。

散点为样本贡献分数：
C_i = 2mean_j||x_i−y_j|| − sum_(j≠i)||x_i−x_j||/(n−1) − Umean(Y,Y)。
C_i的均值等于组能量距离，C_i/n为等权贡献。

共同指标按Normal均值和标准差统一尺度。研究记录按Shannon共同样本规则筛选：熵面板完整、NMD有效且每组每访视至少两条记录。正常参照保留原拟合规则。置换使用筛选后的固定记录集合。
能量距离衡量两组分布在位置、分散程度和形状上的差异。
""",
        "en": """E_U = 2mean||X−Y|| − Umean||X−X′|| − Umean||Y−Y′||.
Y is the Normal reference. Within-group terms use pairs of distinct samples.
Finite-sample estimates can be negative.

Sample contribution scores:
C_i = 2mean_j||x_i−y_j|| − sum_(j≠i)||x_i−x_j||/(n−1) − Umean(Y,Y).
Their mean equals the group energy distance; C_i/n is the equally weighted contribution.

Common markers use Normal means and standard deviations. Study records follow the Shannon common-sample rule: a complete entropy panel, valid NMD and at least two records per group and visit. Normal-reference fitting retains its existing rules. Permutations use the fixed eligible record set.
Energy distance measures differences in location, spread and distribution shape.
""",
    },
    "none": {
        "zh": "本次输出NMD的计算结果和图形。\n",
        "en": "This run provides NMD results and figures.\n",
    },
}

INPUT_GUIDE_EN = """BioEntropy: input naming and format

Load clinical visit tables and animal workbooks together. Experiment types are
detected automatically.

Single-time-point data: experiment.xlsx. Use Group and Sample ID as the first two
columns, followed by numeric markers, with one sample per row.

Longitudinal data: experiment-0W.xlsx, experiment-6W.xlsx, etc. Use the original
patient identifier in USUBJID, stable across visits. Trials 105, 106 and 105and106
use fixed panels and scales. Other experiments use their pooled 0W scale, recorded
in the results.

Normal samples: experiment-Normal.xlsx. Normal.xlsx is accepted in a folder with
one experiment. Multiple Normal files can be attached to an experiment and combined.
Normal ranges: experiment-NormalRange.xlsx, with Feature, Normal_min and Normal_max.
Marker names and units must match exactly; ranges are validated on import.

Select an experiment to set the treatment, comparator and Normal reference.
Drag group names or use Up / Down to set display order. Clear Analyze selected
experiment to skip it. Treatment and comparator settings define comparisons.
The Normal reference determines NMD and energy distance. Unavailable estimates
are recorded as NA with the reason.

Input checks and inclusion rules handle missing or invalid records in the analysis
copy. The input workbooks retain their original contents. English mode translates
the interface, figures, workbook explanations and results viewer. Source marker
headers, group names, patient identifiers and file names retain their original text.

Results:
结果表.xlsx: estimates, sample values, tests, inclusion counts and settings.
图表: figures. 绘图数据: plotting workbooks organized by experiment and figure number.
查看结果.html: results viewer. 本次分析设置.json: saved project. 运行日志.txt: progress log.
"""


def method_guide(method, language="zh"):
    title = METHOD_TITLES[language][method]
    label = "本次方法" if language == "zh" else "Selected method"
    return (
        f"BioEntropy {SOFTWARE_VERSION}\n\n{label}: {title}\n\n"
        f"{COMMON_METHOD_GUIDE[language]}\n## {title}\n\n"
        f"{METHOD_DEFINITIONS[method][language]}"
    )
