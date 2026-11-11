from common_samples import common_records, mark_common
import numpy as np
import pandas as pd
from pathlib import Path
from clinical_inputs import InputSelectionError
from clinical_patient_ids import flat_read
from experiment_config import META, EXTERNAL, NONE
from reference_data import read_reference, clinical_nmd
from entropy_methods import estimate, fixed_bandwidth


def prepare(e):
    frames = []
    columns = None
    arm_by_id = {}
    baseline_ids = set()
    for w, p in zip(e.weeks, e.files):
        h, rr, _ = flat_read(p)
        headers = [str(x).strip() for x in h]
        missing = [k for k in ["Group", "Sample ID", "USUBJID"] if k not in headers]
        if missing:
            raise InputSelectionError(Path(p).name + "：纵向数据缺少列：" + ",".join(missing))
        d = pd.DataFrame([r for _, r in rr], columns=headers)
        d["Excel_row"] = [n for n, _ in rr]
        features = [x for x in headers if x.lower() not in META and x.lower() != "x"]
        if columns is None:
            columns = features
        if set(features) != set(columns):
            raise InputSelectionError(e.name + "：各访视指标或单位不一致。")
        for row in d.itertuples(index=False, name=None):
            key = row[headers.index("USUBJID")]
            group = str(row[headers.index("Group")]).strip()
            if not isinstance(key, str) or not key.strip():
                raise InputSelectionError(Path(p).name + "：USUBJID需为真实患者编号文本。")
            if key in arm_by_id and arm_by_id[key] != group:
                raise InputSelectionError(e.name + "：同一患者跨访视分组不一致。")
            arm_by_id[key] = group
            if w == 0:
                baseline_ids.add(key)
            elif key not in baseline_ids:
                raise InputSelectionError(e.name + "：患者在0W记录中不存在。")
        if d.USUBJID.duplicated().any():
            raise InputSelectionError(Path(p).name + "：同访视患者编号重复。")
        raw = d[columns]
        x = raw.apply(pd.to_numeric, errors="coerce")
        if (raw.notna() & x.isna()).any().any() or np.isinf(x.to_numpy(float)).any():
            raise InputSelectionError(Path(p).name + "：指标含非数值或无穷值。")
        d[columns] = x
        d["Group"] = d.Group.astype(str).str.strip()
        d["Week"] = w
        d["Source_file"] = Path(p).name
        frames.append(d)
    full = pd.concat(frames, ignore_index=True)
    base = full[full.Week.eq(0)][columns]
    features = sorted(c for c in columns if base[c].notna().mean() >= 0.7 and base[c].nunique() > 1)
    if len(features) < 2:
        raise InputSelectionError(e.name + "：基线可用且有变异的共同指标不足2个。")
    center = base[features].mean()
    scale = base[features].std(ddof=1)
    z = (full[features] - center) / scale
    complete = np.isfinite(z).all(1)
    n = int((complete & full.Week.eq(0)).sum())
    if n < 2:
        raise InputSelectionError(e.name + "：基线完整样本不足2个。")
    if e.reference == EXTERNAL:
        ref = read_reference(e.normal_files, features)
    elif e.reference == NONE:
        ref = np.empty((0, len(features)))
    else:
        ref = full.loc[full.Week.eq(0) & full.Group.eq(e.reference), features].to_numpy(float)
    distances, normal = clinical_nmd(
        full[features].to_numpy(float)[:, None, :], list(range(len(features))), ref, e.reference
    )
    return dict(
        full=full,
        features=features,
        center=center,
        scale=scale,
        z=z,
        complete=complete,
        nmd=distances[:, 0],
        normal=normal,
        bandwidth=fixed_bandwidth(n, len(features)),
    )


def analyze(e, out, progress, method="Shannon"):
    from analysis_engine import save_tables, write_json, sha, VERSION

    d = prepare(e)
    df = d["full"]
    paired = method == "Shannon" and d["normal"].get("status") == "可计算"
    original_distance, original_complete = d["nmd"].copy(), d["complete"].copy()
    if paired:
        d["nmd"], d["complete"] = common_records(d["nmd"], d["complete"], df.Group, df.Week)
    rows = []
    obs = []
    for g in e.groups:
        for w in e.weeks:
            present = df.Group.eq(g) & df.Week.eq(w)
            good = present & d["complete"]
            nm = present & np.isfinite(d["nmd"])
            progress(f"{e.name} {g} {w}W：通用固定基线尺度，完整样本{int(good.sum())}")
            rows.append(
                dict(
                    Experiment=e.name,
                    Group=g,
                    Week=w,
                    N_present=int(present.sum()),
                    N_NMD=int(nm.sum()),
                    N_entropy=int(good.sum()),
                    NMD=float(d["nmd"][nm].mean()) if nm.any() else np.nan,
                    p_entropy=len(d["features"]),
                    Bandwidth=d["bandwidth"],
                    **estimate(
                        d["z"].loc[good].to_numpy(float), d["bandwidth"], "clinical", method
                    ),
                )
            )
    for i, row in df.iterrows():
        obs.append(
            dict(
                Experiment=e.name,
                Group=row.Group,
                Week=row.Week,
                USUBJID=row.USUBJID,
                ID_verified=True,
                NMD=d["nmd"][i],
                Included_entropy=bool(d["complete"][i]),
                Source_file=row.Source_file,
                Excel_row=row.Excel_row,
                Missing_entropy_features="; ".join(
                    f for f in d["features"] if not np.isfinite(row[f])
                ),
            )
        )
    if paired:
        mark_common(obs, original_distance, original_complete, d["complete"])
    scales = [
        dict(
            Metric=method,
            Feature=f,
            Center=d["center"][f],
            SD=d["scale"][f],
            Source="本实验0W所有组的共同均值/样本SD",
        )
        for f in d["features"]
    ]
    sources = [
        dict(path=str(Path(p).resolve()), sha256=sha(p)) for p in [*e.files, *e.normal_files]
    ]
    save_tables(out, rows, obs, [], scales, {})
    write_json(
        out / "计算配置与来源.json",
        dict(
            version=VERSION,
            mode="generic_longitudinal",
            features=d["features"],
            normal=d["normal"],
            sources=sources,
            bandwidth=d["bandwidth"],
            shared_shannon_sample=bool(paired),
            scale_source="pooled baseline; observed>=70%; fixed across visits and groups",
            patient_adjustment=False,
        ),
    )
    return dict(
        name=e.name,
        mode="clinical",
        summary=rows,
        observations=obs,
        sources=sources,
        folder=out.name,
        shared_shannon_sample=bool(paired),
        _z=d["z"].to_numpy(float),
        _data=dict(
            values=df[d["features"]].to_numpy(float),
            features=d["features"],
            groups=list(df.Group),
            weeks=list(df.Week),
        ),
    )
