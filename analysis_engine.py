from pathlib import Path
import json, hashlib
from common_samples import common_records, clinical_entropy_panel, mark_common
import numpy as np
import pandas as pd
from clinical_inputs import WEEKS, InputSelectionError
from clinical_patient_ids import load_patient_data, family
from animal_inputs import read_animal, nmd_animal
from entropy_methods import estimate, fixed_bandwidth

VERSION = "V1.0.2"
BASE = Path(__file__).resolve().parent


def plain(v):
    if isinstance(v, np.ndarray):
        return plain(v.tolist())
    if isinstance(v, np.generic):
        return plain(v.item())
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, float) and (not np.isfinite(v)):
        return None
    if isinstance(v, dict):
        return {str(k): plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    return v


def write_json(path, value):
    Path(path).write_text(
        json.dumps(plain(value), ensure_ascii=False, indent=2, allow_nan=False), "utf-8"
    )


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def profiles():
    return json.loads((BASE / "profiles.json").read_text("utf-8"))


def csv(value, path):
    df = value if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
    df.to_csv(path, index=False, encoding="utf-8-sig", na_rep="NA", float_format="%.15g")
    return df


def save_tables(out, summary, observations, exclusions, scales, extras):
    out.mkdir(parents=True, exist_ok=True)
    tables = {
        "组结果": pd.DataFrame(summary),
        "逐样本NMD及纳入": pd.DataFrame(observations),
        "缺失与记录核查": pd.DataFrame(exclusions),
        "固定尺度": pd.DataFrame(scales),
    }
    for name, value in extras.items():
        tables[name] = pd.DataFrame(value)
    for name, df in tables.items():
        csv(df, out / (name + ".csv"))
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter

    book = Workbook()
    book.remove(book.active)
    for name, df in tables.items():
        sheet = book.create_sheet(name[:31])
        sheet.append(list(df.columns) or ["记录"])
        for row_index, row in enumerate(df.itertuples(index=False, name=None), 2):
            vals = [
                None if pd.isna(v) else v.item() if isinstance(v, np.generic) else v for v in row
            ]
            sheet.append(vals)
            for column_index, value in enumerate(vals, 1):
                if isinstance(value, str):
                    sheet.cell(row_index, column_index).data_type = "s"
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for cell in sheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="245674")
            cell.alignment = Alignment(wrap_text=True)
        for j, col in enumerate(df.columns, 1):
            sheet.column_dimensions[get_column_letter(j)].width = min(
                38, max(14, len(str(col)) + 3)
            )
    book.save(out / "结果表.xlsx")


def clinical(bundle, out, progress, options=None):
    options = options or {}
    exp = bundle.experiment
    progress(f"{exp}：读取患者编号及正常参照，核对输入记录")
    d = load_patient_data(exp, profiles()[exp], bundle, options)
    group_names = d["group_names"]
    selected_weeks = options.get("weeks", WEEKS)
    cfg = json.loads((BASE / "clinical_entropy_scale.json").read_text("utf-8"))
    x, jj = clinical_entropy_panel(d, cfg["features"])
    z = (x - np.array(cfg["center"])) / np.array(cfg["scale"])
    good = np.isfinite(z).all(2)
    original_distance, original_complete = d["NMD"].copy(), good.copy()
    paired = (
        options.get("_method", "Shannon") == "Shannon" and d["normal"].get("status") == "可计算"
    )
    if paired:
        d["NMD"], good = common_records(d["NMD"], good, d["groups"][:, None], WEEKS)
    rows, observations, exclusions = ([], [], [])
    for arm, group in enumerate(group_names):
        arm_mask = (d["groups"] == arm) & d["verified"]
        for k, week in enumerate(WEEKS):
            if week not in selected_weeks:
                continue
            keep = arm_mask & good[:, k]
            nm = arm_mask & np.isfinite(d["NMD"][:, k])
            progress(f"{exp} {group} {week}W：计算所选方法，完整人数 {int(keep.sum())}")
            result = estimate(
                z[keep, k], cfg["bandwidth"], "clinical", options.get("_method", "Shannon")
            )
            rows.append(
                dict(
                    Experiment=exp,
                    Group=group,
                    Week=week,
                    N_present=int((arm_mask & d["present"][:, k]).sum()),
                    N_NMD=int(nm.sum()),
                    N_entropy=int(keep.sum()),
                    N_paired_baseline=int((keep & good[:, 0]).sum()),
                    NMD=float(d["NMD"][nm, k].mean()) if nm.any() else np.nan,
                    p_entropy=len(jj),
                    Bandwidth=cfg["bandwidth"],
                    **result,
                )
            )
    for i, k in zip(*np.where(d["present"])):
        row = d["row_info"][i, k]
        missing = [cfg["features"][j] for j in range(len(jj)) if not np.isfinite(x[i, k, j])]
        observations.append(
            dict(
                Experiment=exp,
                Trial=d["study"][i],
                USUBJID="" if not d["verified"][i] else d["ids"][i],
                ID_verified=bool(d["verified"][i]),
                Group=group_names[d["groups"][i]],
                Week=WEEKS[k],
                NMD=d["NMD"][i, k],
                Included_entropy=bool(good[i, k]),
                Source_file=Path(row["file"]).name,
                Excel_row=row["row"],
                Missing_entropy_features="; ".join(missing),
            )
        )
        if missing or not d["verified"][i] or (paired and not good[i, k]):
            exclusions.append(
                dict(
                    USUBJID=d["ids"][i],
                    Group=group_names[d["groups"][i]],
                    Week=WEEKS[k],
                    Reason="编号未核实且化验全缺失"
                    if not d["verified"][i]
                    else ("当次熵面板不完整" if missing else "共同样本纳入条件未满足"),
                    Missing_features="; ".join(missing),
                    Source_file=Path(row["file"]).name,
                    Excel_row=row["row"],
                )
            )
        for j in np.flatnonzero(d["technical"][i, k]):
            exclusions.append(
                dict(
                    USUBJID=d["ids"][i],
                    Group=group_names[d["groups"][i]],
                    Week=WEEKS[k],
                    Reason="负测量值或HbA1c>100置缺失",
                    Feature=d["features"][j],
                    Raw_value=d["Y"][i, k, j],
                    Source_file=Path(row["file"]).name,
                    Excel_row=row["row"],
                )
            )
    if paired:
        mark_common(
            observations,
            original_distance[d["present"]],
            original_complete[d["present"]],
            good[d["present"]],
        )
    scales = [
        dict(
            Metric=options.get("_method", "Shannon"),
            Feature=f,
            Center=cfg["center"][j],
            SD=cfg["scale"][j],
            Source=cfg["scale_id"],
        )
        for j, f in enumerate(cfg["features"])
    ]
    scales += [
        dict(
            Metric="NMD",
            Feature=d["profile"]["nmd_features"][j],
            Center=d["normal"]["mean"][j],
            SD=d["normal"]["sd"][j],
            Source=d["normal"].get("reference", "本次Normal.xlsx"),
        )
        for j in range(len(d["nmd_indices"]))
    ]
    extras = {
        "说明": [
            dict(Item="所选熵方法", Meaning="仅计算所选方法；每组每访视一个估计值"),
            dict(Item="编号", Meaning="核对同一患者、分组和试验来源"),
            dict(Item="数值积分误差", Meaning="Shannon_numeric_SE仅反映计算精度"),
            dict(Item="指标数", Meaning="NMD采用设定面板，两种熵共用固定九指标"),
        ]
    }
    save_tables(out, rows, observations, exclusions, scales, extras)
    meta = dict(
        version=VERSION,
        mode="clinical",
        experiment=exp,
        entropy_scale=cfg,
        normal=d["normal"],
        nmd_features=d["profile"]["nmd_features"],
        sources=d["sources"],
        entropy_integration="Sobol 2**15 nodes x 4 scrambles; seed 20261029"
        if options.get("_method", "Shannon") == "Shannon"
        else "analytic Gaussian kernel overlaps; no quadrature"
        if options.get("_method") == "Renyi2"
        else "not computed",
        shared_shannon_sample=bool(paired),
        missing_entropy="complete cases; no imputation",
        original_NMD_reader_sha256=sha(BASE / "clinical_patient_ids.py"),
        patient_adjustment=False,
        group_order=group_names,
        entropy_estimand="finite equal-weight Gaussian KDE; includes self pairs"
        if options.get("_method") == "Renyi2"
        else "finite equal-weight Gaussian KDE"
        if options.get("_method", "Shannon") == "Shannon"
        else "not computed",
        entropy_uses_NMD=False,
    )
    write_json(out / "计算配置与来源.json", meta)
    return dict(
        name=exp,
        mode="clinical",
        summary=rows,
        observations=observations,
        sources=d["sources"],
        folder=out.name,
        shared_shannon_sample=bool(paired),
        _z=z[d["present"]],
        _data=dict(
            values=x[d["present"]],
            features=cfg["features"],
            groups=[group_names[d["groups"][i]] for i, k in zip(*np.where(d["present"]))],
            weeks=[WEEKS[k] for i, k in zip(*np.where(d["present"]))],
        ),
    )


def animal(path, out, progress, options=None):
    options = options or {}
    from experiment_config import NONE, EXTERNAL
    from reference_data import read_reference

    progress(Path(path).name + "：读取动物数据与正常组")
    before = sha(path)
    d = read_animal(path, options.get("reference"))
    exp = Path(path).stem
    refpaths = options.get("normal_files", [])
    if options.get("reference") == EXTERNAL:
        d["reference_x"] = read_reference(refpaths, d["features"])
    try:
        if options.get("reference") == NONE:
            raise InputSelectionError("未指定正常参照")
        nmd, status, normal = nmd_animal(d)
        reference_ready = True
    except InputSelectionError as exc:
        reference_ready = False
        nmd = np.full(len(d["x"]), np.nan)
        status = np.full(len(d["x"]), "NA：" + str(exc), object)
        normal = dict(
            mean=[np.nan] * len(d["features"]),
            observed_SD_ddof1=[np.nan] * len(d["features"]),
            status="NA：" + str(exc),
        )
        progress("NMD未计算：" + str(exc) + "；两种熵继续计算。")
    paired = options.get("_method", "Shannon") == "Shannon" and reference_ready
    original_distance, original_complete = nmd.copy(), d["complete"].to_numpy().copy()
    if paired:
        nmd, keep = common_records(nmd, d["complete"], d["groups"], 0)
        d["complete"] = pd.Series(keep, index=d["x"].index)
        status = np.where(d["complete"], status, "共同样本纳入条件未满足")
    n = int(d["complete"].sum())
    p = len(d["features"])
    h = fixed_bandwidth(n, p) if options.get("_method", "Shannon") != "none" and n >= 2 else np.nan
    rows, observations, exclusions = ([], [], [])
    for group in d["group_order"]:
        gm = d["groups"].eq(group).to_numpy()
        good = gm & d["complete"].to_numpy()
        nm = gm & np.isfinite(nmd)
        progress(f"{exp} {group}：计算所选方法，完整动物 {int(good.sum())}")
        result = estimate(
            d["z"].loc[good].to_numpy(float), h, "animal", options.get("_method", "Shannon")
        )
        rows.append(
            dict(
                Experiment=exp,
                Group=group,
                Week=0,
                N_present=int(gm.sum()),
                N_NMD=int(nm.sum()),
                N_entropy=int(good.sum()),
                NMD=float(nmd[nm].mean()) if nm.any() else np.nan,
                p_entropy=p,
                Bandwidth=h,
                **result,
            )
        )
    for i, (idx, x) in enumerate(d["x"].iterrows()):
        missing = x.index[x.isna()].tolist()
        observations.append(
            dict(
                Experiment=exp,
                Group=d["groups"].loc[idx],
                Sample_ID=d["ids"].loc[idx],
                NMD=nmd[i],
                NMD_status=status[i],
                Included_entropy=bool(d["complete"].iloc[i]),
                Missing_entropy_features="; ".join(missing),
                Excel_row=int(idx) + 2,
                Source_file=Path(path).name,
            )
        )
        if missing or (paired and not d["complete"].iloc[i]):
            exclusions.append(
                dict(
                    Group=d["groups"].loc[idx],
                    Sample_ID=d["ids"].loc[idx],
                    Missing_features="; ".join(missing),
                    Reason="共同样本纳入条件未满足"
                    if paired
                    else "熵完整病例规则；NMD按部分缺失规则处理",
                )
            )
    if paired:
        mark_common(observations, original_distance, original_complete, d["complete"])
    scales = [
        dict(
            Metric=options.get("_method", "Shannon"),
            Feature=f,
            Center=d["center"][f],
            SD=d["scale"][f],
            Source="本实验所有组的观测均值及样本标准差",
        )
        for f in d["features"]
    ]
    scales += [
        dict(
            Metric="NMD",
            Feature=f,
            Center=normal["mean"][j],
            SD=normal["observed_SD_ddof1"][j],
            Source=d["normal"],
        )
        for j, f in enumerate(d["features"])
    ]
    extras = {
        "未纳入指标": d["excluded_features"],
        "说明": [
            dict(Item="熵", Meaning="基于当组完整样本估计；仅计算所选方法"),
            dict(Item="NMD", Meaning="Normal参照的全协方差马氏距离"),
            dict(Item="数值积分误差", Meaning="Shannon_numeric_SE是计算精度"),
            dict(Item="主面板", Meaning="动物指标可用性规则：非缺失≥70%、有变异；各组相同指标"),
        ],
    }
    save_tables(out, rows, observations, exclusions, scales, extras)
    source = {"path": str(Path(path).resolve()), "sha256": before}
    sources = [source] + [{"path": str(Path(p).resolve()), "sha256": sha(p)} for p in refpaths]
    meta = dict(
        version=VERSION,
        mode="animal",
        experiment=exp,
        normal=normal,
        entropy_features=d["features"],
        entropy_scale=dict(center=d["center"].tolist(), scale=d["scale"].tolist()),
        bandwidth=h,
        bandwidth_n=n,
        sources=sources,
        group_order=d["group_order"],
        entropy_integration="per-kernel Sobol 2**16 nodes x 8 scrambles; seed 20261030"
        if options.get("_method", "Shannon") == "Shannon"
        else "analytic Gaussian kernel overlaps; no quadrature"
        if options.get("_method") == "Renyi2"
        else "not computed",
        shared_shannon_sample=bool(paired),
        missing_entropy="complete cases; no imputation",
        entropy_uses_NMD=False,
        entropy_estimand="finite equal-weight Gaussian KDE; includes self pairs"
        if options.get("_method") == "Renyi2"
        else "finite equal-weight Gaussian KDE"
        if options.get("_method", "Shannon") == "Shannon"
        else "not computed",
        log_transform=False,
        clipping=False,
        winsorization=False,
        contraction=False,
    )
    write_json(out / "计算配置与来源.json", meta)
    return dict(
        name=exp,
        mode="animal",
        summary=rows,
        observations=observations,
        sources=sources,
        folder=out.name,
        shared_shannon_sample=bool(paired),
        _z=d["z"].to_numpy(float),
        _data=dict(
            values=d["x"].to_numpy(float),
            features=d["features"],
            groups=list(d["groups"]),
            weeks=[0] * len(d["x"]),
        ),
    )
