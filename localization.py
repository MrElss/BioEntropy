from contextlib import contextmanager
from contextvars import ContextVar
from html.parser import HTMLParser
import html
import json
import re

LANGUAGE = ContextVar("bioentropy_language", default="zh")
IDENTITIES = ContextVar("bioentropy_identities", default=())

TEXT = {
    "组别": "Group",
    "周": "Week",
    "共同样本清单": "Common sample ledger",
    "个体编号": "Subject ID",
    "共同纳入": "Included in common sample",
    "Excel行号": "Excel row",
    "原因": "Reason",
    "纳入": "Included",
    "缺失指标": "Missing markers",
    "来源文件": "Source file",
    "所选熵面板存在缺失指标": "Incomplete selected entropy panel",
    "NMD不可计算": "NMD unavailable",
    "本组本访视共同有效记录不足2例": "Fewer than two common records in group/visit",
    "V1.0.2": "V1.0.2",
    "BioEntropy 桌面版 · ": "BioEntropy Desktop · ",
    "Shannon熵（默认）": "Shannon entropy (default)",
    "Rényi二阶熵": "Rényi entropy (order 2)",
    "能量距离（正常参照）": "Energy distance (Normal reference)",
    "不计算（仅NMD）": "NMD only",
    "标准检验（9,999次置换）": "Standard test (9,999 permutations)",
    "快速预览（999次置换）": "Quick preview (999 permutations)",
    "关闭显著性检验": "Significance testing off",
    "（未指定）": "(Not specified)",
    "（上传的Normal样本）": "(Uploaded Normal samples)",
    "临床/多时间点": "Clinical / longitudinal",
    "动物/单时间点": "Animal / single time point",
    "计算所选实验": "Analyze selected experiment",
    "文件名": "File",
    "类型": "Type",
    "实验号": "Experiment",
    "实验": "Experiment",
    "周数/说明": "Week / note",
    "记录数": "Records",
    "说明": "Note",
    "计算": "Run",
    "治疗组": "Treatment",
    "比较组/模型组": "Comparator / model",
    "正常参照": "Normal reference",
    "Normal样本": "Normal samples",
    "正常范围": "Normal ranges",
    "参照": "Reference",
    "含真实编号列": "Contains verified patient ID column",
    "未提供USUBJID": "USUBJID not supplied",
    "正常范围仅作核查": "Ranges are for review only",
    "独立正常样本": "Independent Normal samples",
    "请选择数据，程序自动识别实验。": "Select input data. Experiment types are detected automatically.",
    "程序已启动。可同时添加临床分周表和动物分组表；选择实验后可修改组别和显示顺序。": (
        "Ready. Add clinical visit tables and animal workbooks together. "
        "Select an experiment to edit roles and display order."
    ),
    "识别完成。核对组别后即可开始。": "Detection complete. Check the group roles before starting.",
    "输入核查完成": "Input validation complete",
    "核查通过。原始文件保持原样。": "Validation passed. Original files are unchanged.",
    "输入核查通过。": "Input validation passed.",
    "正在分析…": "Analyzing…",
    "分析未完成，请查看日志。": "Analysis did not finish. See the run log.",
    "完成。可点击“查看结果图”或“打开结果”。": "Complete. Use View figures or Open results.",
    "已恢复上次组别与显示顺序；开始分析前会重新核查文件。": "Previous roles and display order restored. Inputs will be checked again before analysis.",
    "正在计算": "Analysis running",
    "请等待当前计算完成。": "Please wait for the current analysis to finish.",
    "关闭将中断计算，是否关闭？": "Closing will interrupt the analysis. Close anyway?",
    "结果目录": "Output folder",
    "文件识别": "Input detection",
    "输入核查": "Input validation",
    "选择实验": "Select experiment",
    "载入设置": "Load project",
    "输入": "Input",
    "请先添加数据。": "Please add input data first.",
    "请选择结果目录。": "Please select an output folder.",
    "请先载入并选中实验。": "Please load and select an experiment first.",
    "选择保存位置（将在其下新建实验结果文件夹）": "Choose a parent folder for the new experiment results",
    "选择实验Excel文件（可多选）": "Select experiment Excel files (multiple files allowed)",
    "选择数据目录（识别本目录及下一层文件夹）": "Choose a data folder (this folder and one level below)",
    "保存分析设置": "Save project",
    "载入分析设置": "Load project",
    "项目设置": "Project settings",
    "命名与填写说明": "Input naming and format",
    "本次方法与核查规则": "Selected method and validation rules",
    "导出模板与填写说明": "Export input templates",
    "选择模板类型": "Template type",
    "临床/参照模板的指标样式": "Marker layout for clinical / reference templates",
    "动物数据": "Animal data",
    "临床分周数据": "Clinical visit data",
    "Normal范围": "Normal ranges",
    "输入填写说明": "Input guide",
    "通用": "Generic",
    "导出": "Export",
    "数据": "Data",
    "选择…": "Choose…",
    "移除所选参照（不删除原文件）": "Remove selected reference (original files are retained)",
    "样本 | ": "Samples | ",
    "范围 | ": "Ranges | ",
    "Hedges’ g值计算": "Hedges’ g calculation",
    "G值数据目录": "G-value data folder",
    "输出目录": "Output folder",
    "开始G值计算": "Calculate G values",
    "G值输入": "G-value input",
    "请选择原G值数据目录。": "Please select the G-value input folder.",
    "单指标、跨模型指标及NMD的G值计算与森林图。": "G values and forest plots for individual endpoints, cross-model endpoints and NMD.",
    "结果预览": "Results preview",
    "此目录没有结果图。": "This folder contains no result figures.",
    "选择已生成的结果目录": "Choose an existing results folder",
    "BioEntropy · 结果图": "BioEntropy · Result figures",
    "本次方法": "Selected method",
    "打开当前图": "Open current figure",
    "打开结果文件夹": "Open results folder",
    "三维交互图：点击右侧“打开当前图”可拖动旋转。": "Interactive 3D figure: choose Open current figure to rotate it.",
    "Shannon熵": "Shannon entropy",
    "能量距离": "Energy distance",
    "仅NMD": "NMD only",
    "NMD（正常参照马氏距离）": "NMD (Normal-reference Mahalanobis distance)",
    "（nats）": " (nats)",
    "基线对齐显示": "Baseline-aligned display",
    "归一化显示": "Normalized display",
    "原始值": "Raw values",
    "基线对齐": "Baseline-aligned",
    "全部散点": "Full scatter distribution",
    " 显示指数（0–100）": " display index (0–100)",
    " 显示指数": " display index",
    "各组减自身0W后共用线性尺度；起点重叠仅用于显示。": "Each group is centered on its own 0W and uses the same linear scale. Baseline alignment is for display only.",
    "按各组估计值的共同范围线性映射。": "Linear mapping uses the common range of group estimates.",
    "原始估计值。": "Raw estimates.",
    " 超出主图的散点用箭头计数，小图保留全部分布。": " Arrows count points outside the main panel; the inset shows the full distribution.",
    "小点：每个样本的NMD。": "Small points: sample NMD values.",
    "小点：样本信息分数Cᵢ；大点/柱高：组Shannon熵。": "Small points: sample information scores Cᵢ; large points / bars: group Shannon entropy.",
    "小点：样本能量距离贡献分数；均值为组估计。": "Small points: sample energy-distance contributions; their mean is the group estimate.",
    "每组/访视一个所选方法的组估计。": "One estimate of the selected metric per group and visit.",
    "两组相对0W变化之差": "between-group difference in change from 0W",
    "当次两组原值之差": "between-group difference in raw values",
    "标记：": "Annotations: ",
    "。探索性置换，Holm校正；†需谨慎核查。": ". Exploratory permutation test, Holm-adjusted; review † cautions.",
    "同一实验的两种指标 · 各自线性显示尺度": "Two metrics from the same experiment · metric-specific linear display scales",
    "同指标组统计量与对应样本分数共用线性变换。完整分布见小图。": "Group estimates and sample scores for each metric share a linear transform. Insets show the full distributions.",
    "各指标采用各自的线性显示尺度。标记仅针对指定两组；临床为两组相对0W变化之差，动物为原值之差；探索性置换、Holm校正。": (
        "Each metric has its own linear display scale. Annotations compare "
        "the selected groups: changes from 0W for clinical data, raw "
        "values for animals; exploratory permutation tests with Holm "
        "adjustment."
    ),
    "每点对应一个组/时间点；用于查看两项指标的联合变化。": "Each point represents a group / visit. The plot displays the joint changes of both metrics.",
    "：NMD 与 ": ": NMD and ",
    " 状态图": " state plot",
    "对照图": "comparison",
    "NMD与": "NMD and ",
    "状态图": "state plot",
    "（显示值）": " (display values)",
    "（原始值）": " (raw values)",
    "未指定": "Not specified",
    "未检验": "Not tested",
    "0W对齐": "0W aligned",
    "NMD（原始值）": "NMD (raw values)",
    "趋势放大版": "Trend-focused view",
    "兼顾散点分布版": "Scatter-distribution view",
    "组统计量视图": "Group-estimate view",
    " · 原始值 →": " · raw values →",
    "完整分布\nNMD原始值": "Full distribution\nNMD (raw values)",
    "：分周 NMD–": ": weekly NMD–",
    " 联合状态图": " joint state plot",
    " · 原始数值 · 各周共用坐标 · 0W不平移、不对齐": " · raw values · shared axes across visits · no 0W shift or alignment",
    "小点：同一患者同次访视的NMD与能量距离贡献分数。": "Small points: NMD and energy-distance contribution from the same patient and visit.",
    "小点：同一患者同次访视的NMD与信息分数。": "Small points: NMD and information score from the same patient and visit.",
    "本图显示所选方法的组统计量。": "This view shows group estimates of the selected metric.",
    "大点：组NMD与所选方法组估计。": "Large points: mean NMD and the selected group estimate. ",
    "箭头为各方向超出主图的数量；角落外的同一点可能被两个方向计数。下方小图完整保留散点。": (
        "Arrows count out-of-window points by direction; a corner point "
        "may count twice. Lower panels retain all scatter points."
    ),
    "每个组点表示当周原始组估计。": "Group points are the raw estimates at each visit.",
    "P值：": "P values: ",
    "，当周两组原始数值之差（含0W）；探索性患者轨迹置换，Holm校正；†见说明。": (
        ", between-group difference in raw values at each visit (including "
        "0W); exploratory patient-trajectory permutation, Holm-adjusted; † "
        "see cautions."
    ),
    "分周联合图 · ": "Weekly joint plot · ",
    "分周原值": "Weekly raw values",
    "常规": "Standard",
    "原安慰剂组": "Original placebo arm",
    " / 样本贡献分数": " / sample contribution score",
    " / 样本信息分数": " / sample information score",
    "本次所选临床实验、两组、全部访视共用线性范围": "Shared linear range across selected trials, groups and visits",
    "组结果": "Group results",
    "样本结果": "Sample results",
    "显著性": "Significance",
    "原值与显示值": "Raw and display values",
    "纳入人数": "Included counts",
    "各图P值": "Figure P values",
    "图表索引": "Figure index",
    "运行设置": "Run settings",
    "输入文件": "Input files",
    "运行": "Run",
    "计算设置": "Calculation settings",
    "统计设置": "Statistical settings",
    "组值": "Groups",
    "散点": "Points",
    "窗口计数": "Window counts",
    "P值": "P values",
    "查看结果": "View results",
    "打开结果表": "Open results table",
    "全部实验": "All experiments",
    "绘图数据": "Plot data",
    "打开图片": "Open image",
    "原始数值比较": "Raw-value comparison",
    "相对0W变化比较": "Change-from-0W comparison",
    "G值": "G values",
    "单指标G值": "Individual endpoint G",
    "模型内指标": "Within-model endpoints",
    "模型内合并": "Within-model pooled",
    "跨模型各模型": "Cross-model estimates",
    "跨模型合并": "Cross-model pooled",
    "NMD汇总": "NMD summary",
    "NMD_G值": "NMD G values",
    "样本NMD": "Sample NMD",
    "输入核查：": "Input validation: ",
    "单时间点分组数据": "Single-time-point grouped data",
    "多时间点/患者数据": "Longitudinal / patient data",
    "：按所选方法计算，并按设定顺序排列图表": ": calculating the selected method and applying the chosen display order",
    "：读取患者编号及正常参照，核对输入记录": ": reading patient IDs and Normal reference; validating input records",
    "：读取动物数据与正常组": ": reading animal data and Normal reference",
    "生成主结果图及绘图数据": "Generating result figures and plot data",
    "整理结果图、绘图数据及结果表": "Organizing figures, plot data and results table",
    "全部完成：": "Complete: ",
    "计算G值并生成森林图": "Calculating G values and generating forest plots",
    "运行未完成；详细原因见错误信息.txt。": "Run incomplete; see 错误信息.txt for details.",
    "请选择空的结果目录。": "Choose an empty output folder.",
    "结果目录应在原始数据目录之外，避免混入输入。": "Choose an output folder outside the input-data folder.",
    "G值结果目录应位于原始G值数据目录之外。": "Choose a G-value output folder outside the G-value input folder.",
    "G值输入文件在运行期间发生变化。": "A G-value input file changed during analysis.",
    "NA：": "NA: ",
    "NA（": "NA (",
    "）：": "): ",
}

PATTERNS = [
    (
        r"(\d+)W：计算样本信息分数（(\d+)个）",
        r"\1W: calculating sample information scores (\2 records)",
    ),
    (
        r"(\d+)W：通用固定基线尺度，完整样本(\d+)",
        r"\1W: fixed baseline scale; complete records: \2",
    ),
    (r"第(\d+)行", r"row \1"),
    (
        r"^已识别 (\d+) 个实验；(\d+) 份实验数据。请选择实验核对组别。$",
        r"Detected \1 experiments and \2 data files. Select an experiment to check group roles.",
    ),
    (
        r"^Normal样本：(\d+)份；正常范围：(\d+)份。正常范围用于核查；NMD采用Normal样本。\n治疗/比较组用于组间比较；正常参照用于NMD；显示顺序只改变图表排列。$",
        (
            "Normal sample files: \\1; range files: \\2. Ranges support "
            "validation.\\nTreatment and comparator define comparisons; the "
            "Normal reference defines NMD. Display order changes presentation "
            "only."
        ),
    ),
    (r"^(\d+)组$", r"\1 groups"),
    (r"，(\d+)组，(\d+)个时间点$", r"; \1 groups; \2 visits"),
    (r"：计算所选方法，完整人数 (\d+)", r": estimating selected metric; complete patients: \1"),
    (r"：计算所选方法，完整动物 (\d+)", r": estimating selected metric; complete animals: \1"),
    (r"：能量距离，纳入(\d+)；", r": energy distance; included: \1; "),
    (
        r" (\d+)W：(\d+)次组别置换；仅检验所选方法",
        r" \1W: \2 group-label permutations; testing selected metric only",
    ),
    (r"^图(\d+)([_ ·])", r"Fig\1\2"),
    (r"^图(\d+)  ", r"Fig\1  "),
    (r"^图(\d+)$", r"Fig\1"),
]

PREFIXES = {
    "开始计算时会自动创建结果目录：": "The output folder will be created when analysis starts: ",
    "已保存设置：": "Project saved: ",
    "已导出：": "Exported: ",
    "结果将保存到：": "Results will be saved to: ",
    "计算完成：": "Analysis complete: ",
    "结果已保存：": "Results saved: ",
    "无法创建结果目录：": "Unable to create output folder: ",
    "：已指定参照文件，请运行前检查输入。": ": reference files assigned; validate inputs before running.",
    "：参照文件": ": reference files",
    "：选择": ": select ",
    "请选择": "Please select ",
    "请至少勾选一个实验。": "Select at least one experiment.",
    "单时间点": "Single time point",
    "多时间点": "Longitudinal",
    "组统计量": "group estimate",
    "原值": "raw values",
    "；": "; ",
    "：": ": ",
    "，": ", ",
    "。": ". ",
    "（": " (",
    "）": ")",
}

from translation_catalog import TEXT as EXTRA_TEXT

TEXT.update(EXTRA_TEXT)


def current_language():
    return LANGUAGE.get()


def set_identities(values):
    return tuple(
        sorted(
            {
                str(value)
                for value in values
                if value is not None and re.search(r"[\u3400-\u9fff]", str(value))
            },
            key=len,
            reverse=True,
        )
    )


@contextmanager
def presentation(language="zh", identities=()):
    if language not in ("zh", "en"):
        raise ValueError("Unsupported interface language: " + str(language))
    a = LANGUAGE.set(language)
    b = IDENTITIES.set(set_identities(identities))
    try:
        yield
    finally:
        IDENTITIES.reset(b)
        LANGUAGE.reset(a)


def translate(value, language=None, identities=None):
    if not isinstance(value, str):
        return value
    protected = IDENTITIES.get() if identities is None else set_identities(identities)
    if (language or current_language()) != "en":
        return value
    if not re.search(r"[\u3400-\u9fff]", value):
        return value
    if value in protected:
        return value
    if value in TEXT:
        return TEXT[value]
    for pattern, replacement in PATTERNS:
        if re.search(pattern, value):
            value = re.sub(pattern, replacement, value)
    substitutions = {}

    def save(match):
        token = f"@@KEEP{len(substitutions)}@@"
        substitutions[token] = match.group(0) if hasattr(match, "group") else match
        return token

    value = re.sub(r'(?i)(?:[a-z]:[\\/]|\\\\)[^\r\n<>"|:]*', save, value)
    value = re.sub(
        r'[^\\/\s:：;；,，<>"\']+\.(?:xlsx|xls|json|csv|html|txt|png|md|exe|ico)\b',
        save,
        value,
        flags=re.I,
    )
    for identity in protected:
        if identity in value:
            value = value.replace(identity, save(identity))
    catalog = {**TEXT, **PREFIXES}
    for source in sorted(catalog, key=len, reverse=True):
        if source in value:
            value = value.replace(source, catalog[source])
    for token, original in substitutions.items():
        value = value.replace(token, original)
    return value


def sheet_title(name):
    return translate(name)[:31]


def translate_reference(value):
    if isinstance(value, str) and "#" in value and current_language() == "en":
        filename, sheet = value.split("#", 1)
        return filename + "#" + sheet_title(sheet)
    return value


def localize_frame(frame):
    frame = frame.copy()
    identity_columns = {
        "Group",
        "Treatment",
        "Comparator",
        "Reference",
        "Normal_reference",
        "Experiment",
        "Subject",
        "USUBJID",
        "Sample_ID",
        "Sample ID",
        "个体编号",
        "实验",
        "组别",
        "来源文件",
        "缺失指标",
        "Trial",
        "Source_file",
        "Source_path",
        "Source_feature",
        "path",
        "Feature",
        "feature",
        "File",
        "Data_file",
        "Figure_file",
    }
    raw_settings = {
        "features",
        "nmd_features",
        "entropy_features",
        "files",
        "normal_files",
        "range_files",
        "groups",
        "order",
    }

    def structured(value, key=""):
        if (
            key
            in {
                "name",
                "experiment",
                "file",
                "Source_file",
                "USUBJID",
                "Sample_ID",
                "Group",
                "Feature",
                "Trial",
                "reference",
                "treatment",
                "comparator",
            }
            or key in raw_settings
        ):
            return value
        if isinstance(value, dict):
            return {field: structured(item, field) for field, item in value.items()}
        if isinstance(value, list):
            return [structured(item) for item in value]
        return translate(value)

    def cell(value, setting=""):
        if setting in raw_settings:
            return value
        if isinstance(value, str) and value.lstrip().startswith(("{", "[")):
            try:
                return json.dumps(structured(json.loads(value)), ensure_ascii=False, allow_nan=True)
            except (ValueError, TypeError):
                pass
        return translate(value)

    for column in frame:
        if str(column) in identity_columns:
            continue
        if frame[column].dtype == object or str(frame[column].dtype).startswith("string"):
            if str(column) == "Value" and "Setting" in frame:
                frame[column] = [
                    cell(value, str(setting))
                    for value, setting in zip(frame[column], frame["Setting"])
                ]
            else:
                frame[column] = frame[column].map(
                    lambda value: (
                        translate_reference(value)
                        if str(column) in {"Data", "P_data", "P_source"}
                        else cell(value)
                    )
                )
    frame.columns = [translate(str(column)) for column in frame.columns]
    return frame


def install_rendering():
    from matplotlib.text import Text
    from matplotlib.figure import Figure

    if getattr(Text.set_text, "_bioentropy_localized", False):
        return
    original = Text.set_text

    def set_text(self, value):
        return original(self, translate(value))

    set_text._bioentropy_localized = True
    Text.set_text = set_text
    original_save = Figure.savefig

    def savefig(self, *args, **kwargs):
        if current_language() == "en":
            for axis in self.axes:
                for artist in axis.texts:
                    value = artist.get_text()
                    if (
                        value.startswith("NMD: ")
                        and "\n" in value
                        and artist.get_transform() is axis.transAxes
                    ):
                        x, y = artist.get_position()
                        if y == 1.02:
                            artist.set_position((x, 1.05))
        return original_save(self, *args, **kwargs)

    Figure.savefig = savefig


class LocalizedHTML(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.parts = []
        self.raw = 0

    def handle_starttag(self, tag, attrs):
        text = self.get_starttag_text()
        if tag in ("script", "style"):
            self.raw += 1
        if tag == "html":
            text = text.replace('lang="zh-CN"', 'lang="en"')
        for key, value in attrs:
            if key in ("title", "alt", "aria-label") and value:
                text = text.replace(
                    key + '="' + html.escape(value, quote=True) + '"',
                    key + '="' + html.escape(translate(value), quote=True) + '"',
                )
        self.parts.append(text)

    def handle_startendtag(self, tag, attrs):
        self.parts.append(self.get_starttag_text())

    def handle_endtag(self, tag):
        self.parts.append("</" + tag + ">")
        if tag in ("script", "style"):
            self.raw = max(0, self.raw - 1)

    def handle_data(self, data):
        self.parts.append(data if self.raw else html.escape(translate(data), quote=False))

    def handle_entityref(self, name):
        self.parts.append("&" + name + ";")

    def handle_charref(self, name):
        self.parts.append("&#" + name + ";")

    def handle_decl(self, decl):
        self.parts.append("<!" + decl + ">")

    def handle_comment(self, data):
        self.parts.append("<!--" + data + "-->")


def localize_html(value):
    if current_language() != "en":
        return value
    parser = LocalizedHTML()
    parser.feed(value)
    return "".join(parser.parts)
