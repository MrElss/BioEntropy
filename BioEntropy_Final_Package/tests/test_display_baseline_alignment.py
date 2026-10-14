"""Display-only contracts: shared anchors, explicit baselines and CSV parity."""
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import matplotlib.pyplot as plt

import bioentropy_display_state as dstate
import bioentropy_figure_data as figure_data
import bioentropy_plots as plots


def _longitudinal_case(metric="NMD"):
    summary = pd.DataFrame({
        "Experiment": ["Clinical"] * 4,
        "Group": ["Treatment", "Treatment", "Placebo", "Placebo"],
        "Week": [0, 6, 0, 6], "n": [2, 2, 2, 2],
        f"{metric}_mean": [12.0, 6.0, 10.0, 8.0],
        "Entropy": [3.5, 2.0, 3.0, 2.6],
        "Entropy_method": ["normal_reference_state_entropy"] * 4,
    }, index=[9, 5, 7, 1])
    samples = pd.DataFrame({
        "Experiment": ["Clinical"] * 8,
        "Group": ["Treatment"] * 4 + ["Placebo"] * 4,
        "Week": [0, 0, 6, 6, 0, 0, 6, 6],
        metric: [10.0, 14.0, 4.0, 8.0, 8.0, 12.0, 6.0, 10.0],
    })
    return summary, samples


@pytest.mark.parametrize("metric", ["NMD_mean", "Entropy"])
def test_different_raw_baselines_align_at_zero_without_mutation(metric):
    summary, _ = _longitudinal_case()
    before = summary.copy(deep=True)
    changes = plots._baseline_change_table(summary, metric)
    baseline = changes[changes["Week"] == 0]
    assert baseline["Baseline_value"].nunique() == 2
    assert baseline["Change_from_baseline"].tolist() == [0.0, 0.0]
    assert baseline["Baseline_reference"].tolist() == ["0W", "0W"]
    np.testing.assert_allclose(changes["Change_from_baseline"],
                               changes[metric] - changes["Baseline_value"])
    assert_frame_equal(summary, before)
    assert changes.index.equals(summary.index)


def test_missing_zero_week_uses_and_labels_each_groups_earliest_valid_week():
    summary, _ = _longitudinal_case()
    summary.loc[summary["Group"] == "Placebo", "Week"] = [6, 12]
    changes = plots._baseline_change_table(summary, "NMD_mean")
    placebo = changes[changes["Group"] == "Placebo"]
    assert placebo["Baseline_week"].tolist() == [6.0, 6.0]
    assert placebo["Baseline_reference"].eq("earliest_available_week").all()
    assert placebo["Change_from_baseline"].tolist() == [0.0, -2.0]
    assert not ((placebo["Week"] == 0).any())


def test_one_timepoint_does_not_invent_longitudinal_change():
    summary, _ = _longitudinal_case()
    one_week = summary[summary["Week"] == 0]
    changes = plots._baseline_change_table(one_week, "NMD_mean")
    assert changes["Change_from_baseline"].isna().all()
    assert changes["Baseline_reference"].eq("unavailable_single_timepoint").all()


def test_baselines_do_not_leak_between_experiments():
    summary, _ = _longitudinal_case()
    second = summary.assign(Experiment="Other", NMD_mean=summary["NMD_mean"] + 50)
    combined = pd.concat([summary, second], ignore_index=True)
    changes = plots._baseline_change_table(combined, "NMD_mean")
    np.testing.assert_allclose(changes["Change_from_baseline"],
                               [0, -6, 0, -2, 0, -6, 0, -2])


def test_no_change_figure_data_when_groups_have_disjoint_single_timepoints(tmp_path):
    summary, samples = _longitudinal_case()
    summary = summary.iloc[[0, 3]]
    results = {"experiments": {"Clinical": {"summary": summary, "sample_level": samples,
                                              "true_normal_available": True}}}
    paths = figure_data.export_figure_data_tables(results, tmp_path)
    assert not any("baseline_change" in p.name for p in paths)


@pytest.mark.parametrize("mode", ["linear", "display"])
def test_scatter_uses_group_summary_anchors_not_individual_extremes(mode):
    summary, samples = _longitudinal_case()
    before = summary.copy(deep=True)
    points = figure_data._collect_points(samples, "Clinical", "NMD",
                                        ["Treatment", "Placebo"], [0, 6], True,
                                        normalize=mode, reference_summary=summary)
    # Group-summary fallback limits are 6 and 12, not individual limits 4 and 14.
    for (group, week), observed in points.items():
        raw = samples.loc[(samples["Group"] == group) & (samples["Week"] == week), "NMD"]
        expected = np.clip((raw.to_numpy() - 6.0) / 6.0, 0.0, 1.0)
        if mode == "display":
            expected = expected ** 3
        np.testing.assert_allclose(observed, expected * 100.0)
    assert_frame_equal(summary, before)


def test_normal_model_anchors_share_equal_group_weights_despite_unequal_sample_counts():
    summary = pd.DataFrame({"Group": ["Normal", "Normal", "Model", "Model"],
                            "Week": [0, 6, 0, 6],
                            "NMD_mean": [1.0, 3.0, 9.0, 11.0]})
    samples = pd.DataFrame({"Experiment": ["Animal"] * 8,
                            "Group": ["Normal"] * 5 + ["Model"] * 3,
                            "Week": [0, 0, 0, 0, 6, 0, 6, 6],
                            "NMD": [1, 1, 1, 1, 3, 9, 11, 11]})
    points = figure_data._collect_points(samples, "Animal", "NMD", ["Normal", "Model"],
                                        [0, 6], True, normalize="linear", reference_summary=summary)
    # Fixed group anchors are Normal=2, Model=10. Individual count weighting is irrelevant.
    np.testing.assert_allclose(points[("Normal", 0)], [0, 0, 0, 0])
    np.testing.assert_allclose(points[("Normal", 6)], [12.5])
    np.testing.assert_allclose(points[("Model", 0)], [87.5])
    np.testing.assert_allclose(points[("Model", 6)], [100, 100])


def test_nonlinearity_does_not_force_mean_of_scatter_to_equal_group_marker(tmp_path):
    summary, samples = _longitudinal_case()
    result = {"experiments": {"Clinical": {"summary": summary, "sample_level": samples,
                                             "true_normal_available": True}}}
    figure_data.export_figure_data_tables(result, tmp_path)
    bars = pd.read_csv(tmp_path / figure_data.FIGURE_DATA_DIRNAME / "Clinical" /
                       "Fig02_NMD_normalized_bars.csv")
    row = bars[(bars["Group"] == "Placebo") & (bars["Week"] == 0)].iloc[0]
    assert row["Display_plotted_0_100_gamma3"] == pytest.approx(100 * (4 / 6) ** 3)
    assert row["Mean_of_plotted_points"] == pytest.approx(50 * ((2 / 6) ** 3 + 1))
    assert row["Bar_equals_mean_of_points"] == "NO"


def test_exported_scatter_matches_plot_overlay(monkeypatch):
    summary, samples = _longitudinal_case()
    monkeypatch.setattr(dstate, "CURRENT_SAMPLE_LEVEL", samples)
    monkeypatch.setattr(dstate, "CURRENT_SUMMARY", summary)
    groups, weeks = ["Treatment", "Placebo"], [0, 6]
    exported = figure_data._collect_points(samples, "Clinical", "NMD", groups, weeks,
                                          True, normalize="display", reference_summary=summary)
    fig, ax = plt.subplots()
    try:
        plotted = plots._overlay_trajectory_points(
            ax, "Clinical", "NMD_mean", weeks, groups,
            {"Treatment": "blue", "Placebo": "red"}, normalize="display")
    finally:
        plt.close(fig)
    np.testing.assert_allclose(plotted, np.concatenate(list(exported.values())))


@pytest.mark.parametrize("method", ["normal_reference_state_entropy", "covariance_entropy"])
def test_entropy_export_has_no_incompatible_individual_dots(tmp_path, method):
    summary, samples = _longitudinal_case()
    summary["Entropy_method"] = method
    result = {"experiments": {"Clinical": {"summary": summary, "sample_level": samples,
                                             "true_normal_available": True}}}
    paths = figure_data.export_figure_data_tables(result, tmp_path)
    names = {p.name for p in paths}
    assert not any("Fig05" in name and "Dots" in name for name in names)
    assert "Fig03_Entropy_raw_bars.csv" in names
    assert "Fig04_Entropy_normalized_bars.csv" in names
    assert not any("Entropy_raw_points" in p or "Entropy_normalized_points" in p or
                   "EntropyDots" in p for p in names)
    mapping = pd.read_csv(tmp_path / figure_data.FIGURE_DATA_DIRNAME / "图表-数据对照表.csv")
    assert not mapping["文件"].str.contains("EntropyDots|Entropy_raw_points|Entropy_normalized_points").any()


@pytest.mark.parametrize("metric", ["NMD", "BRI", "NRBS"])
def test_longitudinal_delta_exports_and_single_timepoint_omission(tmp_path, metric):
    summary, samples = _longitudinal_case(metric)
    before_summary, before_samples = summary.copy(deep=True), samples.copy(deep=True)
    results = {"experiments": {"Clinical": {"summary": summary, "sample_level": samples,
                                              "true_normal_available": metric == "NMD"}}}
    paths = figure_data.export_figure_data_tables(results, tmp_path / "clinical")
    names = {p.name for p in paths}
    for fig, name in (("09", metric), ("10", "Entropy")):
        filename = f"Fig{fig}_{name}_baseline_change_bars.csv"
        assert filename in names
        table = pd.read_csv(tmp_path / "clinical" / figure_data.FIGURE_DATA_DIRNAME /
                            "Clinical" / filename)
        assert table.loc[table["Week"] == 0, "Change_from_baseline"].eq(0).all()
        assert {f"{name}_raw", "Baseline_week", "Baseline_value", "Baseline_reference"}.issubset(table.columns)
    assert_frame_equal(summary, before_summary)
    assert_frame_equal(samples, before_samples)
    results["experiments"]["Clinical"]["summary"] = summary[summary["Week"] == 0]
    results["experiments"]["Clinical"]["sample_level"] = samples[samples["Week"] == 0]
    paths = figure_data.export_figure_data_tables(results, tmp_path / "animal")
    assert not any("baseline_change" in p.name for p in paths)


def test_joint_charts_share_anchors_from_only_valid_joint_rows(tmp_path):
    summary = pd.DataFrame({"Experiment": ["Animal"] * 3, "Week": [0] * 3,
                            "Group": ["Normal", "Model", "Treatment"],
                            "NMD_mean": [2.0, 10.0, 6.0], "Entropy": [np.nan, 3.0, 2.0]})
    samples = pd.DataFrame({"Experiment": ["Animal"] * 4, "Week": [0] * 4,
                            "Group": ["Model", "Model", "Treatment", "Treatment"],
                            "NMD": [8.0, 12.0, 4.0, 8.0]})
    results = {"experiments": {"Animal": {"summary": summary, "sample_level": samples,
                                            "true_normal_available": True}}}
    figure_data.export_figure_data_tables(results, tmp_path)
    folder = tmp_path / figure_data.FIGURE_DATA_DIRNAME / "Animal"
    grouped = pd.read_csv(folder / "Fig05_NMD_Entropy_grouped_bar_NMDDots_points.csv")
    state = pd.read_csv(folder / "Fig08_NMD_Entropy_state_normalized_points.csv")
    grouped_bars = pd.read_csv(folder / "Fig05_NMD_Entropy_grouped_bar_bars.csv")
    assert set(grouped_bars["Group"]) == {"Model", "Treatment"}
    assert sorted(grouped["Value"].tolist()) == [0.0, 12.5, 12.5, 100.0]
    assert sorted(state["Value"].tolist()) == [0.0, 50.0, 50.0, 100.0]
