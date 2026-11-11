from pathlib import Path
from collections import Counter
import copy, hashlib, re
import numpy as np
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
from clinical_inputs import InputSelectionError, WEEKS

GROUPS = ["HTD1801", "安慰剂"]
LAB_NAMES = [
    "γ-谷氨酰转肽酶",
    "丙氨酸氨基转移酶",
    "低密度脂蛋白胆固醇",
    "天门冬氨酸氨基转移酶",
    "总胆固醇",
    "甘油三酯",
    "空腹血浆血糖",
    "糖化血红蛋白",
    "血糖",
    "超敏C反应蛋白",
    "非高密度脂蛋白胆固醇",
    "高密度脂蛋白胆固醇",
]


def numeric(x):
    try:
        y = float(x)
        return y if np.isfinite(y) else np.nan
    except (ValueError, TypeError):
        return np.nan


def flat_read(path):
    if Path(path).suffix.lower() == ".xls":
        import xlrd

        wb = xlrd.open_workbook(str(path), on_demand=True)
        try:
            sh = wb.sheet_by_index(0)

            def values(i):
                return [
                    None
                    if sh.cell_type(i, j) in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK)
                    else sh.cell_value(i, j)
                    for j in range(sh.ncols)
                ]

            h = values(0) if sh.nrows else []
            while h and h[-1] is None:
                h.pop()
            if not h:
                raise InputSelectionError(f"{Path(path).name}：首行没有列名。")
            names = [str(x).strip() for x in h if x is not None and str(x).strip()]
            repeated = [s for s, n in Counter(names).items() if n > 1]
            if repeated:
                raise InputSelectionError(
                    f"{Path(path).name}：存在重复列名：" + "、".join(repeated)
                )
            rows = [(i + 1, values(i)[: len(h)]) for i in range(1, sh.nrows)]
            return (h, [(i, r) for i, r in rows if any(v is not None for v in r)], sh.name)
        finally:
            wb.release_resources()
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        sh = wb.worksheets[0]
        h = list(next(sh.iter_rows(max_row=1, values_only=True)))
        while h and h[-1] is None:
            h.pop()
        if not h:
            raise InputSelectionError(f"{Path(path).name}：首行没有列名。")
        names = [str(x).strip() for x in h if x is not None and str(x).strip()]
        repeated = [s for s, n in Counter(names).items() if n > 1]
        if repeated:
            raise InputSelectionError(f"{Path(path).name}：存在重复列名：" + "、".join(repeated))
        rows = [
            (i, r)
            for i, r in enumerate(sh.iter_rows(min_row=2, max_col=len(h), values_only=True), 2)
            if any(v is not None for v in r)
        ]
        return (h, rows, sh.title)
    finally:
        wb.close()


def family(header):
    h = str(header).strip()
    return next(
        (
            name
            for name in LAB_NAMES
            if h == name
            or h.startswith(name + " ")
            or h.startswith(name + "-")
            or h.startswith(name + "(")
        ),
        None,
    )


def marker_columns(headers, path):
    found = {}
    for j, h in enumerate(headers):
        label = family(h)
        if label:
            if label in found:
                raise InputSelectionError(
                    f"{Path(path).name}：{label}存在多个候选列，无法自动选择。"
                )
            found[label] = j
    missing = [name for name in LAB_NAMES if name not in found]
    if missing:
        raise InputSelectionError(
            f"{Path(path).name}：缺少用于编号及记录核查的化验列：" + "、".join(missing)
        )
    return found


def resolve_trial(rawid, exp, file, row):
    if not isinstance(rawid, str) or not rawid.strip():
        raise InputSelectionError(
            f"{Path(file).name}第{row}行：USUBJID需为原始真实患者编号文本，不能使用表内数字序号。"
        )
    rawid = rawid.strip()
    match = re.search("PCT(105|106)-", rawid)
    if exp == "105and106":
        if not match:
            raise InputSelectionError(
                f"{Path(file).name}第{row}行：合并表的USUBJID需含PCT105-或PCT106-，以区分试验来源。"
            )
        return (rawid, match.group(1))
    if match and match.group(1) != exp:
        raise InputSelectionError(
            f"{Path(file).name}第{row}行：患者编号的试验来源与{exp}文件夹不一致。"
        )
    return (rawid, exp)


def load_patient_data(exp, profile, bundle, options=None):
    options = options or {}
    group_names = options.get("groups", GROUPS)
    visit_weeks = options.get("weeks", WEEKS)
    normal_paths = [
        Path(p)
        for p in options.get("normal_files", [bundle.normal_file] if bundle.normal_file else [])
    ]
    p = copy.deepcopy(profile)
    features = p["features"].copy()
    records = []
    files = list(bundle.visit_files)
    manifest = [
        {"path": str(f.resolve()), "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}
        for f in [*files, *normal_paths]
    ]
    baseline_table = flat_read(files[0])
    base_h = baseline_table[0]
    base_cols = marker_columns(base_h, files[0])
    for name in LAB_NAMES:
        h = str(base_h[base_cols[name]]).strip()
        if not any(family(f) == name for f in features):
            features.append(h)
    baseline = []
    seen = {}
    id_arms = {}
    unresolved = []
    for week, file in zip(visit_weeks, files):
        h, rows, sheet = baseline_table if week == 0 else flat_read(file)
        cols = {str(x).strip(): j for j, x in enumerate(h) if x is not None}
        absent = [f for f in ["Group", "Sample ID", "USUBJID"] if f not in cols]
        if absent:
            raise InputSelectionError(
                f"{file.name}缺少列："
                + "、".join(absent)
                + "。请选择上一轮“已补患者编号”的分周数据文件夹；程序不会用Sample ID代替USUBJID。"
            )
        labs = marker_columns(h, file)
        crp_header = str(h[labs["超敏C反应蛋白"]]).strip()
        if not re.search("\\(\\s*mg/L\\s*\\)", crp_header, re.IGNORECASE):
            raise InputSelectionError(
                f"{file.name}：CRP列需明确标为mg/L；本版CRP阈值按mg/L固定，未自动猜测或换算单位。"
            )
        for f in features:
            if f not in cols and f not in p["features"]:
                name = family(f)
                actual = str(h[labs[name]]).strip() if name in labs else ""
                if name == "超敏C反应蛋白":
                    cols[f] = labs[name]
                elif name == "糖化血红蛋白" and (
                    "(%)" in actual or actual == "糖化血红蛋白-105-ADHBA1C"
                ):
                    cols[f] = labs[name]
        missing = [f for f in features if f not in cols]
        if missing:
            raise InputSelectionError(f"{file.name}缺少固定指标：" + "、".join(missing))
        for row, r in rows:
            arm = str(r[cols["Group"]]).strip()
            if arm not in group_names:
                raise InputSelectionError(f"{file.name}第{row}行：无法识别Group“{arm}”。")
            rawid = r[cols["USUBJID"]]
            verified = rawid is not None and str(rawid).strip() != ""
            values = np.array([numeric(r[cols[f]]) for f in features])
            for f in features:
                raw = r[cols[f]]
                if (
                    raw is not None
                    and str(raw).strip() not in ("", "NA", "N/A", "NaN", "nan")
                    and (not np.isfinite(numeric(raw)))
                ):
                    raise InputSelectionError(
                        f"{file.name}第{row}行：{f}为非数值“{raw}”，请先核查；未自动改成缺失。"
                    )
            if verified:
                key, trial = resolve_trial(rawid, exp, file, row)
            else:
                if np.isfinite(values).any():
                    raise InputSelectionError(
                        f"{file.name}第{row}行：有化验值但USUBJID为空，不能可靠做纵向调整。请补核编号后再运行。"
                    )
                if week == 0:
                    raise InputSelectionError(
                        f"{file.name}第{row}行：基线USUBJID为空，无法确定患者队列。"
                    )
                key = f"UNRESOLVED|{exp}|{week}|{row}"
                trial = ""
                unresolved.append(
                    {
                        "File": file.name,
                        "Sheet": sheet,
                        "Excel_row": row,
                        "Week": week,
                        "Group": arm,
                        "Sample_ID": r[cols["Sample ID"]],
                        "USUBJID": "",
                        "Reason": "12项化验均缺失；编号未确定，保留记录但不参与患者配对",
                    }
                )
            if (key, week) in seen:
                raise InputSelectionError(
                    f"{file.name}第{row}行：患者{key}在{week}W重复，首次出现在第{seen[key, week]}行。"
                )
            seen[key, week] = row
            if key in id_arms and id_arms[key] != arm:
                raise InputSelectionError(
                    f"{file.name}第{row}行：同一USUBJID的Group跨访视不一致：{key}。"
                )
            id_arms[key] = arm
            records.append(
                dict(
                    key=key,
                    arm=arm,
                    week=week,
                    values=values,
                    verified=verified,
                    trial=trial,
                    file=str(file),
                    sheet=sheet,
                    row=row,
                    sample_id=r[cols["Sample ID"]],
                    cells=[f"{sheet}!{get_column_letter(cols[f] + 1)}{row}" for f in features],
                )
            )
            if week == 0:
                baseline.append(key)
    baseline_set = set(baseline)
    extra = [r for r in records if r["verified"] and r["key"] not in baseline_set]
    if extra:
        r = extra[0]
        raise InputSelectionError(
            f"{Path(r['file']).name}第{r['row']}行：USUBJID {r['key']}在0W表中不存在，请核查编号或补充真实基线记录。"
        )
    ids = baseline + [r["key"] for r in records if not r["verified"]]
    index = {s: i for i, s in enumerate(ids)}
    n = len(ids)
    Y = np.full((n, 9, len(features)), np.nan)
    present = np.zeros((n, 9), bool)
    g = np.full(n, -1, int)
    study = np.full(n, "", object)
    verified = np.zeros(n, bool)
    loc = {}
    original_ids = {}
    row_info = {}
    for r in records:
        i = index[r["key"]]
        k = WEEKS.index(r["week"])
        a = group_names.index(r["arm"])
        Y[i, k] = r["values"]
        present[i, k] = True
        g[i] = a
        study[i] = r["trial"]
        verified[i] = r["verified"]
        loc[i, k] = (r["file"], r["cells"])
        original_ids[i] = r["key"]
        row_info[i, k] = r
    if exp == "105and106" and set(study[verified]) != {"105", "106"}:
        raise InputSelectionError(
            "合并数据需同时包含可识别来源的105和106患者。请选择单试验文件夹分析单个试验。"
        )
    technical = np.isfinite(Y) & (Y < 0)
    for j, f in enumerate(features):
        if "糖化血红蛋白" in f:
            technical[:, :, j] |= Y[:, :, j] > 100
    clean = Y.copy()
    clean[technical] = np.nan
    center = np.nanmean(clean[:, 0], axis=0)
    scale = np.nanstd(clean[:, 0], axis=0, ddof=1)
    eh = [features.index(f) for f in p["entropy_features"]]
    nj = [features.index(f) for f in p["nmd_features"]]
    if np.any(~np.isfinite(scale[eh])) or np.any(scale[eh] <= 0):
        raise InputSelectionError("固定熵面板的基线标准差不可用，未计算。")
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, np.nan)
    Z = (clean - center) / scale
    from reference_data import read_reference, clinical_nmd
    from experiment_config import EXTERNAL, NONE

    selected = options.get("reference", EXTERNAL)
    if selected == EXTERNAL:
        A = read_reference(normal_paths, p["normal_features"])
    elif selected == NONE:
        A = np.empty((0, len(nj)))
    else:
        A = clean[g == group_names.index(selected), 0][:, nj]
    N, normal_meta = clinical_nmd(clean, nj, A, selected)
    idx = lambda s: next(j for j, f in enumerate(features) if f.startswith(s))
    rel = Y[:, :, idx("非高密度")] - (Y[:, :, idx("总胆固醇")] - Y[:, :, idx("高密度")])
    relationship = {
        "observed_triplets": int(np.isfinite(rel).sum()),
        "max_absolute_error": float(np.nanmax(np.abs(rel))) if np.isfinite(rel).any() else np.nan,
        "non_hdl_excluded": "derived TC-HDL; retained as QC field",
    }
    p.update(
        longitudinal=True,
        id_status="USUBJID column; baseline, visit uniqueness, group and trial consistency validated",
    )
    for source in manifest:
        if hashlib.sha256(Path(source["path"]).read_bytes()).hexdigest() != source["sha256"]:
            raise InputSelectionError(
                "读取期间文件发生变化，请保存Excel后重新计算：" + source["path"]
            )
    return dict(
        exp=exp,
        profile=p,
        Y=Y,
        Z=Z,
        ids=np.array(ids),
        groups=g,
        study=study,
        verified=verified,
        present=present,
        features=features,
        entropy_indices=eh,
        nmd_indices=nj,
        NMD=N,
        locations=loc,
        original_ids=original_ids,
        row_info=row_info,
        unresolved=unresolved,
        technical=technical,
        center=center,
        scale=scale,
        sources=manifest,
        normal=normal_meta,
        group_names=group_names,
        non_hdl=relationship,
    )
