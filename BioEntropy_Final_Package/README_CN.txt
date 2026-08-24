BioEntropy 使用说明（中文）

BioEntropy 是一个自包含的桌面分析工具，用于对分组的生物医学指标数据
（临床化验 panel 或动物生化数据）做多指标“疾病偏离度”与“系统熵”的探索性
量化。方法与公式的权威说明见 METHODS.md。


一、运行方式（需本机安装 Python 3.11+）
1. 安装依赖：
     cd BioEntropy_Final_Package
     python -m pip install -r requirements_desktop.txt
2. 启动桌面版：
     python BioEntropy_Desktop_App.py
   （Windows 也可双击 Start_Desktop.cmd）
3. 构建单文件 EXE（可选）：
   双击 Clean_And_Rebuild.cmd；构建成功后 EXE 位于 dist\BioEntropy_Desktop.exe。
   若对方运行仍是旧结果，通常是 PyInstaller 旧缓存/旧 spec/旧 dist 未清理，
   先运行 Clean_And_Rebuild.cmd 再重新分发即可。


二、输入数据
- 标准化的 .xlsx 工作簿：每列一个组，每行一个受试者。
- 临床纵向文件命名为“实验号-周数W.xlsx”（如 105-12W.xlsx），软件按 Week 显示。
- 动物实验或其他单时间点文件（如 DB鼠.xlsx）会自动作为 0W 单时间点处理，
  横坐标改用给药组/处理组，而不是把 0W 当作时间轴。
- 健康（Normal）参考：可以是数据内的 Normal/健康组，也可以是单独的 Normal
  工作簿（如 105-Normal.xlsx）。常见美制单位（mg/dL）在中位数明显异常时会
  自动换算为 SI 单位（mmol/L）。


三、核心指标
1. NMD（Normal 参考马氏距离）
   样本到健康参考云的稳健封顶距离。按 Normal 组的稳健均值/尺度转为 z 距离，
   并对单指标极端偏离封顶，降低小样本 Normal 方差过小导致的单指标支配。
   原始 Ledoit-Wolf 协方差距离保留在 NMD_covariance_raw。NMD 越低越接近 Normal。
2. 系统熵（Entropy）
   一组多指标分布的对数行列式（高斯）熵，用 Ledoit-Wolf 收缩估计协方差。
   有 Normal 参考时主 Entropy 使用 Normal-reference state entropy；结果表保留
   Entropy_raw、Entropy_robust_covariance、Entropy_reference_state、Entropy_method
   便于追溯。
3. NRBS（无 Normal 组时）
   以已发表成人参考范围为锚点的绝对负担分数；0 表示所有指标都在参考范围内，
   越高表示偏离越大。
4. BRI（既无 Normal 也无参考范围时）
   相对最早时间点（0W）的基线负担指数，不是绝对健康程度。
   负担指标优先级：NMD > NRBS > BRI。


四、输出
figures 文件夹保留 9 个主图/主文件：疾病偏离度真实值、疾病偏离度归一化值、
熵真实值、熵归一化值、负担 × 熵二维柱形图、可旋转 3D HTML、真实值可点击走势
点图、归一化可点击走势点图、药物机制驱动指标图。结果表以单独的 .xlsx 文件导出。


五、解读提醒
- 每个 Week 的 NMD 和 Entropy 都基于该周原始数据独立计算，不是由上一时间点
  递推得到；12W 是直接用 12W 数据与同一套 Normal 参考比较，而非在 6W 上累加。
- 0W 两组虽都未给药，但随机化后仍可能有基线差异，因此两组熵不一定完全重合。
- 归一化展示（0–100，gamma=3 治疗窗口增强）仅用于图形；正式判断请回到原始
  数值、线性锚定值、置信区间和统计检验。
- 3D HTML 的柱高经过归一化，标签和表格才是原始数值。
- 只有存在真实 Normal 参考时，第一张主图才严格解释为 True-Normal NMD；否则使用
  BRI/HDI 作为基线参考，避免把相对基线指标误称为绝对马氏距离。
