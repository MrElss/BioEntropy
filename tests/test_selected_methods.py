import unittest, tempfile, json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from entropy_methods import estimate
from energy_distance import energy_u, fill_group, labelings
from metric_selection import METHODS, validate_method
from experiment_config import discover, save_project, load_project
from analysis_workflow import run_project


class EnergyTests(unittest.TestCase):
    def test_matches_independent_pair_loops_and_score_mean(self):
        x = np.array([[1.0, 2.0], [2.0, 3.0], [4.0, 6.0]])
        y = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 1.0], [2.0, 2.0]])
        cross = sum((np.linalg.norm(a - b) for a in x for b in y)) / (len(x) * len(y))
        within = lambda a: (
            sum(
                (
                    np.linalg.norm(a[i] - a[j])
                    for i in range(len(a))
                    for j in range(len(a))
                    if i != j
                )
            )
            / (len(a) * (len(a) - 1))
        )
        expected = 2 * cross - within(x) - within(y)
        value, scores = energy_u(x, y, True)
        self.assertAlmostEqual(value, expected, places=12)
        self.assertAlmostEqual(scores.mean(), value, places=12)

    def test_reference_self_comparison(self):
        x = np.array([[0.0], [1.0], [2.0]])
        self.assertAlmostEqual(energy_u(x, x), -2 * cdist(x, x).sum() / (3 * 3 * 2), places=12)
        self.assertLess(energy_u(x, x), 0)

    def test_invariance_and_units(self):
        x = np.arange(12).reshape(6, 2).astype(float)
        y = x[:3] + 0.3
        self.assertAlmostEqual(energy_u(x + 100, y + 100), energy_u(x, y), places=12)
        self.assertAlmostEqual(energy_u(x * 3, y * 3), 3 * energy_u(x, y), places=12)

    def test_median_imputation(self):
        x = np.array([[np.nan, np.nan], [1.0, 2.0], [3.0, np.nan], [5.0, 6.0]])
        a, keep = fill_group(x)
        np.testing.assert_array_equal(keep, [False, True, True, True])
        self.assertEqual(a[1, 1], 4.0)
        self.assertIsNone(fill_group(np.array([[1.0, np.nan], [2.0, np.nan]]))[0])

    def test_permutation_refills_each_group_and_matches_brute_force(self):
        x = np.array([[1.0, np.nan], [2.0, 2.0], [5.0, 4.0], [8.0, 6.0], [10.0, 9.0], [12.0, 13.0]])
        y = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 1.0]])
        a = np.array([[1, 1], [1, 0], [1, 1], [0, 0], [0, 1], [0, 0]], float)
        got = labelings(x, y, a)
        for j in range(2):
            for gi in range(2):
                b, _ = fill_group(x[a[:, j].astype(bool) == (gi == 0)])
                self.assertAlmostEqual(got[j, gi], energy_u(b, y), places=12)
        x = np.nan_to_num(x, nan=5.0)
        got = labelings(x, y, a)
        for j in range(2):
            for gi in range(2):
                self.assertAlmostEqual(
                    got[j, gi], energy_u(x[a[:, j].astype(bool) == (gi == 0)], y), places=12
                )

    def test_method_dispatch(self):
        x = np.array([[0.0, 1.0], [2.0, 3.0], [3.0, 1.0]])
        with patch("entropy_methods.shannon_animal", side_effect=AssertionError("unselected")):
            r = estimate(x, 0.7, "animal", "Renyi2")
            self.assertNotIn("Entropy_Shannon", r)
        with patch("entropy_methods.renyi2", side_effect=AssertionError("unselected")):
            r = estimate(x, 0.7, "animal", "Shannon")
            self.assertNotIn("Entropy_Renyi2", r)
        self.assertEqual(estimate(x, 0.7, "animal", "none"), {})

    def test_invalid_method_selection(self):
        with self.assertRaises(ValueError):
            validate_method("both")
        self.assertIn("Energy", METHODS.values())


class SelectedWorkflowTests(unittest.TestCase):
    def test_nmd_only_accepts_partial_animals_without_complete_kde_records(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            inputs = root / "inputs"
            inputs.mkdir()
            rng = np.random.default_rng(10)
            raw = pd.DataFrame(
                {
                    "Group": ["Normal Control"] * 4 + ["Model Control"] * 4 + ["HTD1801"] * 4,
                    "Sample ID": range(12),
                }
            )
            values = rng.normal(5, 1, (12, 4))
            values[np.arange(12), np.arange(12) % 4] = np.nan
            for j in range(4):
                raw["marker" + str(j)] = values[:, j]
            p = inputs / "partial.xlsx"
            raw.to_excel(p, index=False)
            out = root / "out"
            run_project(discover([p]), out, "none", lambda s: None, 0)
            d = pd.read_excel(out / "结果表.xlsx", sheet_name="组结果")
            self.assertTrue(d.NMD.notna().all())
            self.assertNotIn("Bandwidth", d)

    def test_excel_export_retains_values_and_literal_text(self):
        from table_exports import workbook
        from openpyxl import load_workbook

        with tempfile.TemporaryDirectory() as root:
            p = Path(root) / "check.xlsx"
            frame = pd.DataFrame(
                {
                    "Text": ["=SUM(1,2)", "中文", None],
                    "Number": [1.25, np.nan, -2.0],
                    "Flag": [True, False, True],
                }
            )
            workbook(p, {"Data": frame}, translated=False)
            book = load_workbook(p, data_only=False)
            sheet = book["Data"]
            self.assertEqual(sheet["A2"].value, "=SUM(1,2)")
            self.assertEqual(sheet["A2"].data_type, "s")
            self.assertEqual(sheet["A3"].value, "中文")
            self.assertIsNone(sheet["B3"].value)
            self.assertEqual(sheet["B4"].value, -2.0)
            self.assertEqual(sheet["C3"].value, False)
            book.close()

    def test_all_four_modes_export_only_selected_results(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            input_dir = root / "inputs"
            input_dir.mkdir()
            rng = np.random.default_rng(8)
            raw = pd.DataFrame(
                {
                    "Group": ["Normal Control"] * 4 + ["Model Control"] * 4 + ["HTD1801"] * 4,
                    "Sample ID": range(12),
                    "a": rng.normal(3, 1, 12),
                    "b": rng.normal(4, 1, 12),
                }
            )
            path = input_dir / "fixture.xlsx"
            raw.to_excel(path, index=False)
            original = path.read_bytes()
            for method in ("Shannon", "Renyi2", "Energy", "none"):
                output = root / method
                run_project(discover([path]), output, method, lambda s: None, 0)
                df = pd.read_excel(output / "结果表.xlsx", sheet_name="组结果")
                expected = {
                    "Shannon": "Entropy_Shannon",
                    "Renyi2": "Entropy_Renyi2",
                    "Energy": "Energy_distance",
                }.get(method)
                for field in ("Entropy_Shannon", "Entropy_Renyi2", "Energy_distance"):
                    self.assertEqual(field in df, field == expected)
                index = pd.read_excel(output / "结果表.xlsx", sheet_name="图表索引").to_dict(
                    "records"
                )
                directories = {r["Method"] for r in index}
                self.assertEqual(directories, {"NMD"} | ({method} if method != "none" else set()))
                self.assertTrue(all((r["Method"] in directories for r in index)))
                obs = pd.read_excel(output / "结果表.xlsx", sheet_name="样本结果")
                self.assertEqual("Shannon_information_score" in obs, method == "Shannon")
                self.assertEqual("Energy_contribution_score" in obs, method == "Energy")
                tests = pd.read_excel(output / "结果表.xlsx", sheet_name="显著性")
                self.assertEqual(set(tests.Metric), {"NMD"} | ({expected} if expected else set()))
                self.assertEqual(path.read_bytes(), original)

    def test_generic_longitudinal_energy_and_saved_method(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            inputs = root / "inputs"
            inputs.mkdir()
            rng = np.random.default_rng(9)
            paths = []
            for week in (0, 6):
                d = pd.DataFrame(
                    {
                        "Group": ["Normal Control"] * 4 + ["Model Control"] * 4 + ["HTD1801"] * 4,
                        "Sample ID": range(12),
                        "USUBJID": ["P" + str(i) for i in range(12)],
                        "a": rng.normal(5, 1, 12),
                        "b": rng.normal(5, 1, 12),
                    }
                )
                p = inputs / f"fixture-{week}W.xlsx"
                d.to_excel(p, index=False)
                paths.append(p)
            experiments = discover(paths)
            project = root / "settings.json"
            save_project(project, experiments, "Energy", permutations=0)
            experiments, config = load_project(project)
            output = root / "output"
            run_project(experiments, output, config["entropy_display"], lambda s: None, 0)
            df = pd.read_excel(output / "结果表.xlsx", sheet_name="组结果")
            self.assertTrue(df.Energy_distance.notna().all())
            display = pd.read_excel(output / "结果表.xlsx", sheet_name="原值与显示值")
            for _, part in display.groupby("Metric"):
                self.assertEqual(part[part.Week.eq(0)].Display.nunique(), 1)

    def test_energy_reference_requirement(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            inputs = root / "inputs"
            inputs.mkdir()
            p = inputs / "fixture.xlsx"
            pd.DataFrame(
                {
                    "Group": ["A"] * 3 + ["B"] * 3,
                    "Sample ID": range(6),
                    "a": [1, 2, 3, 4, 5, 6],
                    "b": [2, 1, 3, 5, 4, 6],
                }
            ).to_excel(p, index=False)
            out = root / "out"
            run_project(discover([p]), out, "Energy", lambda s: None, 0)
            d = pd.read_excel(out / "结果表.xlsx", sheet_name="组结果")
            self.assertTrue(d.Energy_distance.isna().all())
            self.assertTrue(d.Energy_status.str.startswith("NA").all())


if __name__ == "__main__":
    unittest.main()
