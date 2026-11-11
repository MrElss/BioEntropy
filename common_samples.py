import numpy as np
import pandas as pd


def clinical_entropy_panel(data, features):
    from clinical_patient_ids import family

    indices = [
        next(i for i, f in enumerate(data["features"]) if family(f) == marker)
        for marker in features
    ]
    values = data["Y"][:, :, indices].copy()
    values[data["technical"][:, :, indices]] = np.nan
    values[~data["verified"]] = np.nan
    return values, indices


def mark_common(observations, original_distance, complete, keep):
    for row, distance, full, included in zip(
        observations,
        np.asarray(original_distance).ravel(),
        np.asarray(complete).ravel(),
        np.asarray(keep).ravel(),
        strict=True,
    ):
        reason = (
            "Included"
            if included
            else "Incomplete entropy panel"
            if not full
            else "NMD unavailable"
            if not np.isfinite(distance)
            else "Fewer than two common records in group/visit"
        )
        row.update(
            Included_common=bool(included), Included_NMD=bool(included), Selection_reason=reason
        )


def common_records(distance, complete, group, week):
    distance = np.asarray(distance, dtype=float)
    eligible = (np.isfinite(distance) & np.asarray(complete, dtype=bool)).ravel()
    cells = pd.DataFrame(
        {
            "group": np.broadcast_to(group, distance.shape).ravel(),
            "week": np.broadcast_to(week, distance.shape).ravel(),
            "eligible": eligible,
        }
    )
    count = cells.groupby(["group", "week"])["eligible"].transform("sum")
    keep = (eligible & (count.to_numpy() >= 2)).reshape(distance.shape)
    return np.where(keep, distance, np.nan), keep
