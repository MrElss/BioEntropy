import unittest, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from entropy_methods import renyi2, shannon_animal, estimate
from animal_inputs import read_animal, nmd_animal
from analysis_engine import profiles
from clinical_inputs import InputBundle, WEEKS, InputSelectionError
from clinical_patient_ids import load_patient_data
from plots import aligned_table
from gvalue_recompute import hedges_g


class CalculationTests(unittest.TestCase):
    def test_identical_kernels_are_one_gaussian(self):
        x = np.full((5, 3), 2.0)
        h = 0.7
        self.assertAlmostEqual(renyi2(x, h)[0], 1.5 * np.log(4 * np.pi * h * h), places=12)
        self.assertAlmostEqual(
            shannon_animal(x, h, power=8, reps=2)[0],
            1.5 * np.log(2 * np.pi * np.e * h * h),
            places=12,
        )

    def test_renyi_includes_self_terms(self):
        x = np.array([[0.0], [2.0], [4.0]])
        h = 0.6
        overlap = sum((np.exp(-float((a - b)[0] ** 2) / (4 * h * h)) for a in x for b in x)) / 9
        self.assertAlmostEqual(
            renyi2(x, h)[0], 0.5 * np.log(4 * np.pi * h * h) - np.log(overlap), places=12
        )

    def test_translation_and_joint_scale(self):
        x = np.random.default_rng(5).normal(size=(5, 2))
        h = 0.8
        r = renyi2(x, h)[0]
        self.assertAlmostEqual(r, renyi2(x + 19, h)[0], places=12)
        self.assertAlmostEqual(renyi2(x * 3, h * 3)[0] - r, 2 * np.log(3), places=12)

    def test_entropy_ordering(self):
        x = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 0.0]])
        hs, se = shannon_animal(x, 0.8, power=10, reps=4)
        self.assertGreaterEqual(hs + 5 * se, renyi2(x, 0.8)[0])

    def test_small_group_status(self):
        a = estimate(np.ones((1, 2)), 0.8, "animal")
        self.assertTrue(np.isnan(a["Entropy_Shannon"]) and np.isnan(a["Entropy_Renyi2"]))
        self.assertIn("NA", a["Entropy_status"])

    def test_alignment_is_linear_reversible_and_common(self):
        results = [
            dict(
                name="x",
                mode="clinical",
                summary=[
                    dict(Group=g, Week=w, NMD=v, Entropy_Shannon=v + 2, Entropy_Renyi2=v + 1)
                    for g, arr in [("A", [9, 4, 3]), ("B", [7, 6, 2])]
                    for w, v in zip([0, 6, 12], arr)
                ],
            )
        ]
        table, bounds = aligned_table(results)
        for m, part in table.groupby("Metric"):
            self.assertEqual(part[part.Week == 0].Display.nunique(), 1)
            np.testing.assert_allclose(
                part.Display * bounds[m]["span"] / 100 + bounds[m]["delta_min"] + part.Baseline_raw,
                part.Raw,
            )

    def test_g_matches_manual_hedges(self):
        x = np.array([5, 6, 7, 8, 9.0])
        y = np.array([1, 2, 3, 4.0])
        n = len(x)
        m = len(y)
        sp = np.sqrt(((n - 1) * x.var(ddof=1) + (m - 1) * y.var(ddof=1)) / (n + m - 2))
        g = (1 - 3 / (4 * (n + m - 2) - 1)) * (x.mean() - y.mean()) / sp
        self.assertAlmostEqual(hedges_g(x, y).g, g, places=12)
        self.assertIsNone(hedges_g([2], [3, 4]))

    def test_animal_missing_and_duplicate_rules(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "animal.xlsx"
            df = pd.DataFrame(
                {
                    "Group": ["Normal Control"] * 4 + ["Model Control"] * 4,
                    "Sample ID": list(range(8)),
                    "a": [1, 2, 3, 4, 5, 6, 7, 8],
                    "b": [2, 4, 3, 8, 9, 11, 10, np.nan],
                }
            )
            df.to_excel(path, index=False)
            d = read_animal(path)
            dist, status, meta = nmd_animal(d)
            self.assertFalse(d["complete"].iloc[-1])
            self.assertTrue(np.isfinite(dist[-1]))
            self.assertIn("部分缺失", status[-1])
            df.loc[1, "Sample ID"] = 0
            df.to_excel(path, index=False)
            with self.assertRaisesRegex(InputSelectionError, "重复"):
                read_animal(path)


class PatientTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.profile = profiles()["106"]
        self.files = []
        for week in WEEKS:
            rows = []
            for i in range(8):
                row = {
                    "Group": "HTD1801" if i < 4 else "安慰剂",
                    "Sample ID": i + 1,
                    "USUBJID": f"PCT106-S{i:03}",
                }
                row.update(
                    {
                        f: float(2 + i * 0.2 + j * 0.01 + week * 0.001)
                        for j, f in enumerate(self.profile["features"])
                    }
                )
                rows.append(row)
            df = pd.DataFrame(rows)
            if week > 0:
                df = df.iloc[::-1]
            path = self.root / f"106-{week}W.xlsx"
            df.to_excel(path, index=False)
            self.files.append(path)
        normal = pd.DataFrame(
            {
                f: np.arange(6) / 10 + 1 + j * 0.03
                for j, f in enumerate(self.profile["normal_features"])
            }
        )
        self.normal = self.root / "Normal.xlsx"
        normal.to_excel(self.normal, index=False)
        self.bundle = InputBundle("106", self.root, self.normal, visit_files=tuple(self.files))

    def tearDown(self):
        self.tmp.cleanup()

    def test_true_ids_survive_reordered_rows(self):
        d = load_patient_data("106", self.profile, self.bundle)
        np.testing.assert_allclose(d["Y"][:, 1, 0] - d["Y"][:, 0, 0], 0.006)
        self.assertEqual(len(d["ids"]), 8)

    def test_duplicate_visit_id_stops(self):
        df = pd.read_excel(self.files[1])
        df.loc[1, "USUBJID"] = df.loc[0, "USUBJID"]
        df.to_excel(self.files[1], index=False)
        with self.assertRaisesRegex(InputSelectionError, "重复"):
            load_patient_data("106", self.profile, self.bundle)

    def test_patient_identifier_requirement(self):
        df = pd.read_excel(self.files[1])
        df.loc[0, "USUBJID"] = None
        df.to_excel(self.files[1], index=False)
        with self.assertRaisesRegex(InputSelectionError, "USUBJID为空"):
            load_patient_data("106", self.profile, self.bundle)

    def test_id_arm_conflict_stops(self):
        df = pd.read_excel(self.files[1])
        df.loc[0, "Group"] = "HTD1801"
        df.to_excel(self.files[1], index=False)
        with self.assertRaisesRegex(InputSelectionError, "不一致"):
            load_patient_data("106", self.profile, self.bundle)


if __name__ == "__main__":
    unittest.main()
