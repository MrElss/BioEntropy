from pathlib import Path
import json
import numpy as np
from common_samples import common_records, clinical_entropy_panel, mark_common
from analysis_engine import save_tables, write_json, sha, VERSION, profiles
from clinical_inputs import WEEKS
from clinical_patient_ids import load_patient_data
from animal_inputs import read_animal, nmd_animal
from reference_data import read_reference
from experiment_config import EXTERNAL, NONE
from energy_distance import energy_u, fill_group, prepare_reference


def attach(
    name,
    mode,
    values,
    features,
    reference,
    groups,
    weeks,
    observations,
    nmd,
    group_order,
    selected_weeks,
    sources,
    out,
    normal_meta,
    roles,
    progress,
    entropy_complete,
):
    values = np.asarray(values, float)
    ref, mu, sd, status = prepare_reference(reference)
    raw_z = (values - mu) / sd if ref is not None else np.full_like(values, np.nan)
    groups, weeks = (np.asarray(groups), np.asarray(weeks))
    original_nmd = np.asarray(nmd, float).copy()
    nmd, common = common_records(nmd, entropy_complete, groups, weeks)
    mark_common(observations, original_nmd, entropy_complete, common)
    raw_z[~common] = np.nan
    z = np.full_like(values, np.nan)
    rows, exclusions = ([], [])
    for index, row in enumerate(observations):
        row.update(
            NMD=nmd[index],
            Included_entropy=False,
            Energy_contribution_score=np.nan,
            Energy_weighted_contribution=np.nan,
            Energy_partial_imputed=False,
        )
        if not common[index]:
            row["NMD_status"] = row["Selection_reason"]
    for group in group_order:
        for week in selected_weeks:
            present = (groups == group) & (weeks == week)
            indices = np.flatnonzero(present)
            eligible = np.flatnonzero(present & common)
            filled, keep = fill_group(raw_z[eligible])
            good = eligible[keep] if filled is not None and ref is not None else np.array([], int)
            nvalid = present & np.isfinite(nmd)
            value = np.nan
            group_status = (
                status if ref is None else "NA：可用记录不足2个或组内至少一项指标全部缺失"
            )
            if ref is not None and filled is not None:
                value, scores = energy_u(filled, ref, scores=True)
                z[good] = filled
                group_status = "OK"
                for index, score in zip(good, scores):
                    observations[index].update(
                        Included_entropy=True,
                        Energy_contribution_score=float(score),
                        Energy_weighted_contribution=float(score / len(good)),
                        Energy_partial_imputed=bool(not np.isfinite(values[index]).all()),
                    )
            for index in indices:
                if index not in good or not np.isfinite(values[index]).all():
                    exclusions.append(
                        dict(
                            Group=group,
                            Week=int(week),
                            Excel_row=observations[index].get("Excel_row"),
                            Sample_ID=observations[index].get("Sample_ID"),
                            USUBJID=observations[index].get("USUBJID"),
                            Reason=observations[index]["Selection_reason"]
                            if not common[index]
                            else group_status
                            if index not in good
                            else "能量距离：部分缺失按本组本访视指标中位数填补；原表不变",
                            Missing_features="; ".join(
                                f
                                for j, f in enumerate(features)
                                if not np.isfinite(values[index, j])
                            ),
                        )
                    )
            progress(f"{name} {group} {week}W：能量距离，纳入{len(good)}；{group_status}")
            rows.append(
                dict(
                    Experiment=name,
                    Group=group,
                    Week=int(week),
                    N_present=int(present.sum()),
                    N_NMD=int(nvalid.sum()),
                    N_entropy=len(good),
                    p_entropy=len(features),
                    NMD=float(nmd[nvalid].mean()) if nvalid.any() else np.nan,
                    Energy_distance=value,
                    Energy_status=group_status,
                )
            )
    scales = [
        dict(
            Metric="Energy_distance",
            Feature=f,
            Center=float(mu[j]) if mu is not None else np.nan,
            SD=float(sd[j]) if sd is not None else np.nan,
            Source="指定Normal样本均值/原始观测样本SD",
        )
        for j, f in enumerate(features)
    ]
    extras = {
        "说明": [
            dict(Item="能量距离", Meaning="整体分布与Normal的差异；直接输出U估计，保留负值"),
            dict(Item="散点", Meaning="组能量距离的样本贡献分数；均值等于组估计"),
            dict(
                Item="缺失",
                Meaning="与Shannon使用相同共同记录：熵面板完整、NMD有效且每组每访视至少2例；正常参照保留原拟合规则",
            ),
            dict(Item="患者编号", Meaning="身份核查和轨迹置换"),
        ]
    }
    save_tables(out, rows, observations, exclusions, scales, extras)
    meta = dict(
        version=VERSION,
        selected_method="Energy",
        mode=mode,
        sources=sources,
        normal=normal_meta,
        energy_features=features,
        nmd_features=normal_meta.get("features", features),
        energy_scale=dict(center=mu, scale=sd, status=status),
        energy_estimator="2*mean(group,Normal)-Umean(group,group)-Umean(Normal,Normal); no square root",
        negative_values="retained; including same reference sample",
        energy_uses_NMD=False,
        shared_analysis_sample=True,
        sample_selection="Complete Shannon panel and valid NMD; at least two records per group/visit",
        missing_energy="Study records selected by the common Shannon/NMD eligibility rule; Normal reference retains its fitting rules",
        patient_adjustment=False,
        log_transform=False,
        clipping=False,
        winsorization=False,
        contraction=False,
    )
    write_json(Path(out) / "计算配置与来源.json", meta)
    return dict(
        name=name,
        mode=mode,
        selected_method="Energy",
        summary=rows,
        observations=observations,
        sources=sources,
        folder=Path(out).name,
        roles=roles,
        shared_analysis_sample=True,
        _z=z,
        _energy_raw=raw_z,
        _energy_reference=ref,
        _data=dict(
            values=values, features=features, groups=list(groups), weeks=list(map(int, weeks))
        ),
    )


def clinical_energy(bundle, out, progress, options):
    d = load_patient_data(bundle.experiment, profiles()[bundle.experiment], bundle, options)
    clean = d["Y"].copy()
    clean[d["technical"]] = np.nan
    cfg = json.loads(Path(__file__).with_name("clinical_entropy_scale.json").read_text("utf-8"))
    entropy_values, _ = clinical_entropy_panel(d, cfg["features"])
    entropy_complete = np.isfinite(entropy_values).all(2)
    indices = d["nmd_indices"]
    features = [d["features"][j] for j in indices]
    selected = options.get("reference", EXTERNAL)
    if selected == EXTERNAL:
        reference = read_reference(
            options.get("normal_files") or [bundle.normal_file], d["profile"]["normal_features"]
        )
    elif selected == NONE:
        reference = np.empty((0, len(indices)))
    else:
        reference = clean[d["groups"] == d["group_names"].index(selected), 0][:, indices]
    present_indices = list(zip(*np.where(d["present"])))
    values = np.array(
        [
            clean[i, k, indices] if d["verified"][i] else np.full(len(indices), np.nan)
            for i, k in present_indices
        ]
    )
    observations = []
    for i, k in present_indices:
        row = d["row_info"][i, k]
        observations.append(
            dict(
                Experiment=bundle.experiment,
                Trial=d["study"][i],
                USUBJID=d["ids"][i] if d["verified"][i] else "",
                ID_verified=bool(d["verified"][i]),
                Group=d["group_names"][d["groups"][i]],
                Week=WEEKS[k],
                NMD=d["NMD"][i, k],
                Source_file=Path(row["file"]).name,
                Excel_row=row["row"],
                Missing_entropy_features="; ".join(
                    f
                    for j, f in enumerate(cfg["features"])
                    if not np.isfinite(entropy_values[i, k, j])
                ),
            )
        )
    return attach(
        bundle.experiment,
        "clinical",
        values,
        features,
        reference,
        [r["Group"] for r in observations],
        [r["Week"] for r in observations],
        observations,
        np.array([r["NMD"] for r in observations]),
        d["group_names"],
        options.get("weeks", WEEKS),
        d["sources"],
        out,
        d["normal"],
        options,
        progress,
        entropy_complete[d["present"]],
    )


def animal_energy(path, out, progress, options):
    from clinical_inputs import InputSelectionError

    d = read_animal(path, options.get("reference"))
    if options.get("reference") == EXTERNAL:
        d["reference_x"] = read_reference(options["normal_files"], d["features"])
    reference = d.get("reference_x", d["x"].loc[d["groups"].eq(d["normal"])].to_numpy(float))
    try:
        if options.get("reference") == NONE:
            raise InputSelectionError("未指定正常参照")
        nmd, status, normal = nmd_animal(d)
    except InputSelectionError as exc:
        nmd = np.full(len(d["x"]), np.nan)
        status = ["NA：" + str(exc)] * len(nmd)
        normal = dict(status="NA：" + str(exc))
    observations = [
        dict(
            Experiment=Path(path).stem,
            Group=d["groups"].loc[idx],
            Sample_ID=d["ids"].loc[idx],
            NMD=nmd[i],
            NMD_status=status[i],
            Excel_row=int(idx) + 2,
            Source_file=Path(path).name,
            Missing_entropy_features="; ".join(
                f for f in d["features"] if not np.isfinite(d["x"].loc[idx, f])
            ),
        )
        for i, idx in enumerate(d["x"].index)
    ]
    sources = [
        dict(path=str(Path(p).resolve()), sha256=sha(p))
        for p in [path, *options.get("normal_files", [])]
    ]
    return attach(
        Path(path).stem,
        "animal",
        d["x"].to_numpy(float),
        d["features"],
        reference,
        list(d["groups"]),
        [0] * len(nmd),
        observations,
        nmd,
        d["group_order"],
        [0],
        sources,
        out,
        normal,
        options,
        progress,
        d["complete"].to_numpy(),
    )


def generic_energy(e, out, progress):
    from generic_longitudinal import prepare

    d = prepare(e)
    df = d["full"]
    features = d["features"]
    if e.reference == EXTERNAL:
        reference = read_reference(e.normal_files, features)
    elif e.reference == NONE:
        reference = np.empty((0, len(features)))
    else:
        reference = df.loc[df.Week.eq(0) & df.Group.eq(e.reference), features].to_numpy(float)
    obs = [
        dict(
            Experiment=e.name,
            Group=row.Group,
            Week=row.Week,
            USUBJID=row.USUBJID,
            ID_verified=True,
            NMD=d["nmd"][i],
            Source_file=row.Source_file,
            Excel_row=row.Excel_row,
            Missing_entropy_features="; ".join(f for f in features if not np.isfinite(row[f])),
        )
        for i, row in df.iterrows()
    ]
    sources = [
        dict(path=str(Path(p).resolve()), sha256=sha(p)) for p in [*e.files, *e.normal_files]
    ]
    return attach(
        e.name,
        "clinical",
        df[features].to_numpy(float),
        features,
        reference,
        list(df.Group),
        list(df.Week),
        obs,
        d["nmd"],
        e.groups,
        e.weeks,
        sources,
        out,
        d["normal"],
        e.record(),
        progress,
        d["complete"].to_numpy(),
    )
