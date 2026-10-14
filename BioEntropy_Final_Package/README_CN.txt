BioEntropy V1.0

BioEntropy 用于分组生物医学指标的探索性综合分析，支持动物实验与临床多时间点数据。

一、启动
Windows 发行包：双击 BioEntropy_V1_0.exe。
源码：在 BioEntropy_Final_Package 中运行 Start_Desktop.cmd，需要 Python 3.11 或更新版本。
也可安装 requirements_desktop.txt 后运行 python BioEntropy_Desktop_App.py。
界面显示 BioEntropy V1.0。结果目录前缀为 BioEntropy_Results_V1_0，
本地设置保存在用户目录 .bioentropy_v1_0_desktop.json。

二、输入
每行一个样本：Group | Sample ID | 指标1 | 指标2 ...。
数值型 ID 在同一文件内唯一。多时间点文件按 实验号-周数W.xlsx 命名。
Normal.xlsx 或 实验号-Normal.xlsx 提供正常样本；正常范围文件使用 NormalRange。
请使用随附输入模板，并核对物种、指标单位、分组角色和正常参照来源。

三、计算与输出
主 NMD 使用 Ledoit-Wolf 全协方差马氏距离，待测样本不作 ±8 z 值封顶。
NMD_diagonal_capped 保留 ±8 封顶的对角距离，用作敏感性对照。
异常检测保留为质量控制标记，不按异常投票移动或排除样本。
共同预处理、截尾、缺失值填补和协方差正则化仍按方法说明执行。
Entropy 优先报告有效 NMD 的二阶矩状态分数，不满足条件时报告协方差熵。
Entropy_covariance 单独保留群体协方差熵。距离衍生分数与 NMD 不构成独立证据。
输出包括结果表、图表、方法敏感性表及同级 *_audit 运行记录。
独立森林图模块使用自己的指标筛选、±3 截尾及效应量合并流程，不能混用两套方法描述。

四、复现与说明
requirements_lock.txt 记录 Python 3.12 的固定分析依赖；requirements_build_lock.txt 记录构建依赖。
输入、参数、依赖和字体可能影响复现，运行记录包含输入哈希、计算版本和方法设置。
方法详情见 METHODS.md、METRIC_DICTIONARY.md、ANALYSIS_GUIDE.md，以及中英文补充方法。
公开发行包仅附模板和合成测试，不含研究原始数据。

五、自检
BioEntropy_V1_0.exe --self-test 自检报告.json
源码测试：python -m pytest tests/ -q
