import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from analysis_engine import animal
from animal_inputs import read_animal, nmd_animal
from energy_analysis import animal_energy
from energy_distance import prepare_reference, energy_u
from experiment_config import discover
from comparison_statistics import compare


class SharedMethodTests(unittest.TestCase):
    def test_energy_uses_same_identifiers_distances_and_permutation_units(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rng = np.random.default_rng(34)
            values = rng.normal(5, 1, (18, 3))
            values[[7, 13], 2] = np.nan
            frame = pd.DataFrame(values, columns=["a", "b", "c"])
            frame.insert(0, "Sample ID", list(range(1, 7)) * 3)
            frame.insert(
                0, "Group", ["Normal Control"] * 6 + ["Model Control"] * 6 + ["HTD1801"] * 6
            )
            path = root / "animals.xlsx"
            frame.to_excel(path, index=False)
            original = path.read_bytes()
            options = discover([path])[0].record()
            with patch("analysis_engine.estimate", return_value={"Entropy_Shannon": 1.0}):
                sh = animal(path, root / "sh", lambda x: None, dict(options, _method="Shannon"))
            en = animal_energy(path, root / "en", lambda x: None, options)
            shobs, enobs = pd.DataFrame(sh["observations"]), pd.DataFrame(en["observations"])
            for field in ["NMD", "Included_entropy", "Included_common", "Selection_reason"]:
                pd.testing.assert_series_equal(shobs[field], enobs[field])
            np.testing.assert_allclose(
                [x["NMD"] for x in sh["summary"]], [x["NMD"] for x in en["summary"]], rtol=0, atol=0
            )
            self.assertEqual([x["N_entropy"] for x in en["summary"]], [6, 5, 5])
            raw = read_animal(path)
            reference, mu, sd, _ = prepare_reference(
                raw["x"].loc[raw["groups"].eq(raw["normal"])].to_numpy()
            )
            for row in en["summary"]:
                use = enobs.Group.eq(row["Group"]) & enobs.Included_common
                expected = energy_u((raw["x"].loc[use].to_numpy() - mu) / sd, reference)
                self.assertAlmostEqual(row["Energy_distance"], expected, places=12)
                self.assertAlmostEqual(
                    enobs.loc[use, "Energy_contribution_score"].mean(), expected, places=12
                )
            en["roles"] = {key: options[key] for key in ["treatment", "comparator", "reference"]}
            en["selected_method"] = "Energy"
            with patch(
                "comparison_statistics.assignments",
                wraps=__import__("comparison_statistics").assignments,
            ) as assign:
                compare(en, 999, lambda x: None)
            self.assertEqual(len(assign.call_args.args[0]), 10)
            self.assertEqual(path.read_bytes(), original)

    def test_blank_excel_rows_do_not_change_identity_matching(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            frame = pd.DataFrame(
                {
                    "Group": ["Normal Control"] * 4 + [None] + ["Model Control"] * 4,
                    "Sample ID": [1, 2, 3, 4, None, 1, 2, 3, 4],
                    "a": [1, 2, 4, 3, None, 5, 6, 7, 8],
                    "b": [2, 5, 1, 3, None, 4, np.nan, 8, 9],
                }
            )
            path = root / "gap.xlsx"
            frame.to_excel(path, index=False)
            with patch("analysis_engine.estimate", return_value={}):
                sh = animal(path, root / "sh", lambda x: None, {"_method": "Shannon"})
            en = animal_energy(path, root / "en", lambda x: None, {})
            pd.testing.assert_series_equal(
                pd.DataFrame(sh["observations"]).NMD, pd.DataFrame(en["observations"]).NMD
            )
            excluded = [x for x in en["observations"] if not x["Included_common"]]
            self.assertEqual(excluded[0]["Sample_ID"], 2)
            self.assertEqual(excluded[0]["Excel_row"], 8)
            self.assertEqual(excluded[0]["Missing_entropy_features"], "b")

    def test_single_complete_animal_is_excluded_in_both_modes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            frame = pd.DataFrame(
                {
                    "Group": ["Normal Control"] * 6 + ["Model Control"] * 3,
                    "Sample ID": range(9),
                    "a": range(1, 10),
                    "b": [2, 5, 1, 3, 7, 4, 8, np.nan, np.nan],
                }
            )
            path = root / "singleton.xlsx"
            frame.to_excel(path, index=False)
            with patch("analysis_engine.estimate", return_value={}):
                sh = animal(path, root / "sh", lambda x: None, {"_method": "Shannon"})
            en = animal_energy(path, root / "en", lambda x: None, {})
            self.assertEqual([r["N_NMD"] for r in sh["summary"]], [6, 0])
            self.assertEqual([r["N_entropy"] for r in en["summary"]], [6, 0])
            self.assertEqual(
                en["observations"][6]["Selection_reason"],
                "Fewer than two common records in group/visit",
            )


if __name__ == "__main__":
    unittest.main()
