# BioEntropy V1.0 指标字典

本表定义 BioEntropy V1.0 的主要输出字段。具体方法由方法标识列和运行记录说明。

| 字段 | 本版定义 | 解释边界 |
|---|---|---|
| NMD / NMD_sq | Ledoit-Wolf 全 Normal 协方差马氏距离 / 距离平方 | 主距离不封顶；依赖参考估计和面板 |
| NMD_covariance_raw | 主 NMD 的兼容别名 | raw 不表示未预处理单位 |
| NMD_method | ledoitwolf_full_covariance | 当前默认主算法 |
| NMD_z_cap | 空值 | 主距离不使用 z 封顶 |
| NMD_diagonal_capped | 原 Normal 均值/SD 下的 ±8 封顶标准化距离 | 单独敏感性分数，不能把两种差异只归因于协方差 |
| NMD_p_chi2 | 全协方差主距离平方的 chi-square 尾面积 | 探索性近似，不是精确校准的检验 |
| NMD_diagonal_p_chi2 | 原对角版本参考尾面积 | 仅用于对照 |
| NMD_p_smallsample_F / NMD_cov_p_smallsample_F | 原 predictive F 缩放，m>p+1 时输出 | 正则化与预处理使精确 F 校准不成立 |
| NMD_feature_count | 实际非缺失参考指标数 | 全缺失仍可能为零距离，不能据此解释为健康 |
| Entropy_covariance | 基线预处理坐标中的 Ledoit-Wolf 高斯协方差熵 | 取消额外个体收缩；不等于异常值稳健估计 |
| Entropy_raw / Entropy_robust_covariance | Entropy_covariance 的别名 | raw 已经过共同预处理；robust 是保留的旧字段名称 |
| Entropy_reference_state | 0.5*log(mean(NMD²))+0.5*log(2πe) | 主全协方差距离二阶矩分数，直接依赖 NMD |
| Entropy / Entropy_method | 有>=2有效NMD且分数有限时使用状态分数，否则协方差熵 | 先核对方法；状态分数与 NMD 的关联并非独立证据 |
| Entropy_n_adjusted_subjects | 0 | 不因异常标记移动/排除样本 |
| Entropy_n_high_confidence_outliers | 0，旧收缩人数别名 | 不代表没有 QC 强标记；不要把旧列当作新 QC 计数 |
| Entropy_n_strong_qc_flags | >=4 方法标记人数 | QC 描述；不是经校准的异常概率 |
| NRBS / BRI / HDI | 参考范围或基线负担 | NRBS 需核对物种、单位、人群；BRI 是相对基线负担 |
| Calculation_version | V1.0 | 组时间点结果的方法版本 |
| 森林图的熵效应量 | 独立模块逐个体 log(distance) 的 Hedges' g | 与主组水平二阶矩/协方差熵定义不同，本次未改该模块 |

`NMD方法敏感性` 给出两种距离的组均值和 bootstrap 区间，`治疗比较` 保存两种组间对照。
新增字段和版本写入结果，audit 中保存 method_policy。旧版结果不与本版混合解释。

图形数据中的 `Baseline_week`、`Baseline_value` 和 `Change_from_baseline` 分别记录
每组实际使用的基线周、该周原始指标值及当前值减基线值。优先使用0W，缺失时使用最早有效周；
不足两个有效时间点不定义变化轨迹。这些仅供主图09/10展示，不改写原始指标或增加统计检验。
归一化个体散点和组曲线共用组汇总值拟合的锚点，熵图仅显示组水平指标。
