import unittest, tempfile, copy
from pathlib import Path
import numpy as np
import pandas as pd
from experiment_config import discover, save_project, load_project, EXTERNAL, NONE
from clinical_inputs import InputSelectionError
from animal_inputs import read_animal, nmd_animal
from entropy_methods import renyi2, fixed_bandwidth
from supplementary import review_flags, range_bounds
from plots import data_for_plot, display_table
from reference_data import clinical_nmd
from generic_longitudinal import prepare


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.raw = pd.DataFrame(
            {
                "group": ["Normal Control"] * 5 + ["Model Control"] * 5 + ["HTD1801"] * 5,
                "sample ID": list(range(15)),
                "a": np.arange(15) + 1.0,
                "b": np.random.default_rng(44).normal(10, 2, 15),
            }
        )
        self.path = self.root / "动物.xlsx"
        self.raw.to_excel(self.path, index=False)

    def tearDown(self):
        self.tmp.cleanup()

    def test_auto_detect_case_insensitive_animal_and_roles(self):
        e = discover([self.path])[0]
        self.assertEqual(e.kind, "animal")
        self.assertEqual(e.reference, "Normal Control")
        self.assertEqual(e.comparator, "Model Control")
        self.assertEqual(e.treatment, "HTD1801")

    def test_selected_file_scope(self):
        self.raw.to_excel(self.root / "未选择.xlsx", index=False)
        self.assertEqual([e.name for e in discover([self.path])], ["动物"])

    def test_normal_range_validation(self):
        pd.DataFrame(
            {"Feature": ["a", "b"], "Normal_min": [0, 0], "Normal_max": [10, 20]}
        ).to_excel(self.root / "动物-NormalRange.xlsx", index=False)
        e = discover([self.root])[0]
        self.assertEqual(len(e.range_files), 1)
        self.assertEqual(e.normal_files, [])
        self.assertEqual(e.reference, "Normal Control")

    def test_external_reference_can_be_added_to_nonstandard_groups(self):
        raw = self.raw.copy()
        raw["group"] = raw.group.str.replace("Normal Control", "参照A")
        raw.to_excel(self.path, index=False)
        pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [3.0, 2.0, 4.0]}).to_excel(
            self.root / "动物-Normal.xlsx", index=False
        )
        e = discover([self.root])[0]
        self.assertEqual(e.reference, EXTERNAL)

    def test_ambiguous_reference(self):
        self.raw.to_excel(self.root / "第二实验.xlsx", index=False)
        self.raw.to_excel(self.root / "Normal.xlsx", index=False)
        with self.assertRaisesRegex(InputSelectionError, "无法唯一匹配"):
            discover([self.root])

    def test_project_restores_roles_order_and_exclusion(self):
        e = discover([self.path])[0]
        e.order.reverse()
        e.treatment = "Model Control"
        e.comparator = "HTD1801"
        e.enabled = False
        path = self.root / "project.json"
        save_project(path, [e])
        loaded, _ = load_project(path)
        self.assertEqual(e.record(), loaded[0].record())

    def test_same_treatment_comparator_rejected(self):
        e = discover([self.path])[0]
        e.treatment = e.comparator
        with self.assertRaisesRegex(InputSelectionError, "不能相同"):
            e.validate()

    def test_order_must_be_a_permutation(self):
        e = discover([self.path])[0]
        e.order[-1] = e.order[0]
        with self.assertRaisesRegex(InputSelectionError, "每个组一次"):
            e.validate()

    def test_reference_changes_nmd_not_entropy_input(self):
        a = read_animal(self.path, "Normal Control")
        b = read_animal(self.path, "Model Control")
        self.assertFalse(np.allclose(nmd_animal(a)[0], nmd_animal(b)[0]))
        pd.testing.assert_frame_equal(a["z"], b["z"])
        h = fixed_bandwidth(len(a["z"]), 2)
        self.assertEqual(renyi2(a["z"].to_numpy(), h), renyi2(b["z"].to_numpy(), h))

    def test_animal_marker_selection(self):
        self.raw["USUBJID"] = ["S" + str(i) for i in range(15)]
        self.raw.to_excel(self.path, index=False)
        self.assertEqual(read_animal(self.path)["features"], ["a", "b"])

    def test_missing_normal_reference(self):
        n, meta = clinical_nmd(np.ones((3, 1, 2)), [0, 1], np.empty((0, 2)), NONE)
        self.assertTrue(np.isnan(n).all())
        self.assertIn("NA", meta["status"])

    def test_range_conflict_rejected(self):
        paths = []
        for i, hi in enumerate([5, 7]):
            p = self.root / f"r{i}.xlsx"
            pd.DataFrame({"Feature": ["a"], "Normal_min": [0], "Normal_max": [hi]}).to_excel(
                p, index=False
            )
            paths.append(p)
        with self.assertRaisesRegex(InputSelectionError, "冲突"):
            range_bounds(paths)

    def test_review_flag_measurements(self):
        x = np.array([[1.0, 2], [2.0, 3], [3.0, 4], [4.0, 5], [100.0, 500.0]])
        before = x.copy()
        flags, _ = review_flags(dict(values=x, groups=["g"] * 5, weeks=[0] * 5))
        np.testing.assert_array_equal(x, before)
        self.assertGreater(flags.Flag_count.iloc[-1], 0)

    def test_reordered_rows_only_change_display_order(self):
        rows = [
            dict(
                Experiment="test",
                Group=g,
                Week=0,
                NMD=v,
                Entropy_Shannon=v + 3,
                Entropy_Renyi2=v + 2,
                N_NMD=5,
                N_entropy=5,
            )
            for g, v in [("A", 2), ("B", 4), ("C", 3)]
        ]
        a = dict(name="test", mode="animal", summary=rows)
        b = copy.deepcopy(a)
        b["summary"].reverse()
        da = data_for_plot(a, "NMD", True, display_table([a], pd.DataFrame()))
        db = data_for_plot(b, "NMD", True, display_table([b], pd.DataFrame()))
        self.assertEqual(db.Group.tolist(), ["C", "B", "A"])
        pd.testing.assert_series_equal(
            da.set_index("Group").Display.sort_index(), db.set_index("Group").Display.sort_index()
        )

    def test_legacy_xls_values_match_xlsx(self):
        import xlwt
        from clinical_patient_ids import flat_read

        path = self.root / "旧格式.xls"
        book = xlwt.Workbook()
        sheet = book.add_sheet("数据")
        for j, v in enumerate(self.raw.columns):
            sheet.write(0, j, v)
        for i, row in enumerate(self.raw.itertuples(index=False, name=None), 1):
            for j, v in enumerate(row):
                sheet.write(i, j, v)
        book.save(str(path))
        e = discover([path])[0]
        self.assertEqual(e.kind, "animal")
        a = read_animal(path)
        b = read_animal(self.path)
        np.testing.assert_allclose(nmd_animal(a)[0], nmd_animal(b)[0], rtol=1e-14)
        self.assertEqual(len(flat_read(path)[1]), len(self.raw))

    def test_generic_longitudinal_ids_and_baseline_scale(self):
        for w in [0, 6]:
            d = self.raw.rename(columns={"group": "Group", "sample ID": "Sample ID"}).copy()
            d["USUBJID"] = ["T-" + str(i) for i in range(15)]
            d[["a", "b"]] += w * 0.01
            if w:
                d = d.iloc[::-1]
            d.to_excel(self.root / f"新试验-{w}W.xlsx", index=False)
        e = discover([self.root / "新试验-0W.xlsx", self.root / "新试验-6W.xlsx"])[0]
        d = prepare(e)
        self.assertEqual(e.kind, "clinical")
        self.assertEqual(len(d["nmd"]), 30)
        self.assertTrue(np.isfinite(d["nmd"]).all())
        pd.testing.assert_series_equal(d["center"], self.raw[["a", "b"]].mean(), check_names=False)


if __name__ == "__main__":
    unittest.main()
