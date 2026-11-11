import numpy as np
import pandas as pd
from metric_selection import score_field


def affine(displays, name, metric):
    a = displays[(displays.Experiment.astype(str) == str(name)) & displays.Metric.eq(metric)]
    bases = dict(a[["Group", "Baseline_raw"]].drop_duplicates().itertuples(index=False, name=None))
    finite = a.dropna(subset=["Change", "Display"]).drop_duplicates("Change")
    if len(finite) >= 2:
        u, v = (finite.iloc[0], finite.iloc[-1])
        slope = float((v.Display - u.Display) / (v.Change - u.Change))
        offset = float(u.Display - slope * u.Change)
    else:
        slope, offset = (0.0, 50.0)
    return (bases, slope, offset)


def coordinates(result, metric, displays):
    groups = pd.DataFrame(result["summary"]).copy()
    points = pd.DataFrame(result["observations"]).copy()
    if "Week" not in points:
        points["Week"] = 0
    for axis, key, field in [("X", metric, score_field(metric)), ("Y", "NMD", "NMD")]:
        bases, slope, offset = affine(displays, result["name"], key)
        groups[axis] = (groups[key] - groups.Group.map(bases)) * slope + offset
        if score_field(key) is not None:
            points[axis] = (
                (points[field] - points.Group.map(bases)) * slope + offset
                if field in points
                else np.nan
            )
        else:
            points[axis] = np.nan
    points["Joint_plot_included"] = np.isfinite(points[["X", "Y"]]).all(axis=1)
    return (groups, points)
