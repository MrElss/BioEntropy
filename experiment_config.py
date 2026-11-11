from pathlib import Path
from dataclasses import dataclass, field, asdict
import json, re
from clinical_inputs import InputSelectionError
from clinical_patient_ids import flat_read

EXTERNAL = "（上传的Normal样本）"
NONE = "（无）"
VERSION = "V1.0.2"
VISIT = re.compile("^(.*?)\\s*[-_ ]\\s*(\\d+)\\s*(?:W|周)$", re.I)
META = {
    "group",
    "sample id",
    "sample_id",
    "usubjid",
    "patient id",
    "patient_id",
    "trial",
    "study",
    "week",
}


def is_normal(group):
    s = "".join(str(group).lower().split())
    return s in (
        "normal",
        "normalcontrol",
        "healthy",
        "healthycontrol",
        "sham",
        "blank",
        "naive",
        "正常",
        "正常组",
        "正常对照组",
        "健康组",
        "健康对照组",
        "空白组",
        "假手术组",
        "control",
    )


@dataclass
class Experiment:
    name: str
    kind: str
    files: list
    weeks: list
    groups: list
    normal_files: list = field(default_factory=list)
    range_files: list = field(default_factory=list)
    treatment: str = NONE
    comparator: str = NONE
    reference: str = NONE
    order: list = field(default_factory=list)
    enabled: bool = True
    rows: list = field(default_factory=list)

    def validate(self):
        if not self.enabled:
            return
        if not self.name or self.name in (".", "..") or re.search('[\\\\/:*?"<>|]', self.name):
            raise InputSelectionError("实验名不可作为安全的结果目录名。")
        if (
            self.kind not in ("animal", "clinical")
            or not self.files
            or len(self.files) != len(self.weeks)
            or (len(set(self.weeks)) != len(self.weeks))
        ):
            raise InputSelectionError(self.name + "：实验类型或访视文件配置无效。")
        if not self.groups or len(set(self.groups)) != len(self.groups):
            raise InputSelectionError(self.name + "：分组配置无效。")
        if set(self.order) != set(self.groups) or len(self.order) != len(self.groups):
            raise InputSelectionError(self.name + "：显示顺序应包含每个组一次。")
        for label, value in [("治疗组", self.treatment), ("比较组", self.comparator)]:
            if value not in [NONE, *self.groups]:
                raise InputSelectionError(self.name + "：" + label + "不存在。")
        if self.treatment == self.comparator and self.treatment != NONE:
            raise InputSelectionError(self.name + "：治疗组和比较组不能相同。")
        if self.reference not in [NONE, EXTERNAL, *self.groups]:
            raise InputSelectionError(self.name + "：正常参照选择无效。")
        if self.reference == EXTERNAL and (not self.normal_files):
            raise InputSelectionError(self.name + "：选择了外部正常参照，但未添加Normal样本文件。")
        if self.kind == "clinical" and 0 not in self.weeks:
            raise InputSelectionError(
                self.name + "：纵向数据需包含0W文件，用于核查真实患者编号和固定尺度。"
            )
        for p in [*self.files, *self.normal_files, *self.range_files]:
            if not Path(p).is_file():
                raise InputSelectionError("文件不存在：" + p)

    def record(self):
        return asdict(self)


def candidates(inputs):
    paths = []
    for value in inputs:
        p = Path(value).expanduser().resolve()
        if not p.exists():
            raise InputSelectionError("找不到输入：" + str(p))
        found = sorted(p.glob("*.xls*")) + sorted(p.glob("*/*.xls*")) if p.is_dir() else [p]
        for f in found:
            if (
                f.suffix.lower() in (".xlsx", ".xls")
                and (not f.name.startswith("~$"))
                and (f not in paths)
            ):
                paths.append(f)
    return paths


def classify(path):
    h, rows, _ = flat_read(path)
    names = [str(x).strip() if x is not None else "" for x in h]
    first = {str(r[0]).strip().lower() for _, r in rows if r}
    lower = {x.lower() for x in names}
    range_name = bool(re.search("normal[ _-]*range|正常范围|正常区间", path.stem, re.I))
    isrange = (
        range_name or {"normal_min", "normal_max"} <= lower or {"normal_min", "normal_max"} <= first
    )
    if isrange:
        return ("range", h, rows, None, None)
    m = VISIT.match(path.stem)
    if m:
        return ("clinical", h, rows, m.group(1), int(m.group(2)))
    if re.search(
        "(?:^|[-_ ])(?:normal|healthy|正常|健康)(?:样本|数据)?(?:$|[-_ ]\\d+$)", path.stem, re.I
    ):
        return ("normal", h, rows, None, None)
    if "group" in lower and ("sample id" in lower or "sample_id" in lower):
        return ("animal", h, rows, path.stem, 0)
    raise InputSelectionError(
        path.name
        + "：无法识别数据类型。数据表需含Group、Sample ID；分周文件命名为“实验名-0W.xlsx”；正常样本命名为“实验名-Normal.xlsx”。"
    )


def discover(inputs):
    items = candidates(inputs)
    experiments = {}
    auxiliary = []
    overview = []
    if not items:
        raise InputSelectionError("所选位置未发现.xlsx或.xls文件。")
    for path in items:
        kind, h, rows, name, week = classify(path)
        if kind in ("normal", "range"):
            auxiliary.append((path, kind, h, rows))
            continue
        key = (str(path.parent), name)
        if key not in experiments:
            experiments[key] = Experiment(name, kind, [], [], [])
        e = experiments[key]
        if e.kind != kind:
            raise InputSelectionError(name + "：同一实验混入单时间点和分周命名，请核对。")
        if week in e.weeks:
            raise InputSelectionError(name + "：同一时间点有多个文件，请只选择一份。")
        labels = [str(x).strip() if x is not None else "" for x in h]
        if "group" not in [x.lower() for x in labels]:
            raise InputSelectionError(path.name + "：缺少Group列。")
        gi = [x.lower() for x in labels].index("group")
        gs = []
        for num, row in rows:
            g = str(row[gi]).strip() if row[gi] is not None else ""
            if not g:
                raise InputSelectionError(f"{path.name}第{num}行：组名为空。")
            if g not in gs:
                gs.append(g)
            if g not in e.groups:
                e.groups.append(g)
        e.files.append(str(path))
        e.weeks.append(week)
        e.rows.append(
            dict(
                file=str(path),
                kind=kind,
                experiment=name,
                week=week,
                count=len(rows),
                note="含真实编号列"
                if "USUBJID" in labels
                else "未提供USUBJID"
                if kind == "clinical"
                else f"{len(gs)}组",
            )
        )
    if not experiments:
        raise InputSelectionError("只发现正常样本或正常范围文件，尚未选择实验数据。")
    names = [e.name for e in experiments.values()]
    if len(set(names)) != len(names):
        raise InputSelectionError("不同目录有同名实验，请分次运行，避免合并不同数据版本。")
    for path, kind, h, rows in auxiliary:
        prefix = re.sub(
            "[-_ ]*(normal[ _-]*range|(?:normal|healthy|正常|健康)(?:样本|数据)?|正常范围|正常区间)(?:[-_ ]\\d+)?$",
            "",
            path.stem,
            flags=re.I,
        )
        matches = [e for (folder, name), e in experiments.items() if name == prefix]
        if not matches and prefix == "":
            matches = [e for (folder, name), e in experiments.items() if folder == str(path.parent)]
        if not matches and prefix == "" and (len(experiments) == 1):
            matches = list(experiments.values())
        if len(matches) != 1:
            raise InputSelectionError(
                path.name
                + "：无法唯一匹配所属实验。请命名为“实验名-Normal.xlsx”或先只载入实验数据，再在该实验设置中指定参照文件。"
            )
        e = matches[0]
        getattr(e, "normal_files" if kind == "normal" else "range_files").append(str(path))
        e.rows.append(
            dict(
                file=str(path),
                kind=kind,
                experiment=e.name,
                week="参照",
                count=len(rows),
                note="正常范围仅作核查" if kind == "range" else "独立正常样本",
            )
        )
    for e in experiments.values():
        pairs = sorted(zip(e.weeks, e.files))
        e.weeks = [w for w, _ in pairs]
        e.files = [p for _, p in pairs]
        healthy = [g for g in e.groups if is_normal(g)]
        e.groups = (
            [*healthy, *[g for g in e.groups if g not in healthy]]
            if e.kind == "animal"
            else ["HTD1801", "安慰剂"]
            if set(e.groups) == {"HTD1801", "安慰剂"}
            else e.groups
        )
        e.reference = EXTERNAL if e.normal_files else healthy[0] if len(healthy) == 1 else NONE
        e.comparator = next(
            (g for g in e.groups if re.search("model|modle|placebo|安慰剂|模型", g, re.I)), NONE
        )
        e.treatment = next(
            (g for g in e.groups if re.search("HTD|1801|治疗", g, re.I) and g != e.comparator), NONE
        )
        e.order = e.groups.copy()
    return list(experiments.values())


def save_project(path, experiments, display="Shannon", output="", permutations=9999):
    data = dict(
        schema=1,
        version=VERSION,
        entropy_display=display,
        output=output,
        permutations=permutations,
        experiments=[e.record() for e in experiments],
    )
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), "utf-8")


def load_project(path):
    data = json.loads(Path(path).read_text("utf-8"))
    if data.get("schema") != 1:
        raise InputSelectionError("项目文件格式不支持。")
    result = [Experiment(**x) for x in data["experiments"]]
    if len({e.name for e in result}) != len(result):
        raise InputSelectionError("项目中有重复实验名。")
    for e in result:
        e.validate()
    return (result, data)
