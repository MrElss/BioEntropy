import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import csv
import numpy as np
import pandas as pd
from gvalue_recompute import compute_maha_model, write_part4_csv, two_sided_p
from normal_reference import SOFTWARE_VERSION


class GReferenceTests(unittest.TestCase):
    def test_reference_exports_and_complete_nmd_g_pipeline(self):
        self.assertEqual(SOFTWARE_VERSION, "V1.0.2")
        rng = np.random.default_rng(317)
        rows = []
        for group, shift in [("Control", 0.0), ("Model", 2.0), ("PM", 1.4), ("HTD1801", 0.8)]:
            for i, values in enumerate(10.0 + shift + rng.normal(0.0, 0.2, (8, 3))):
                rows.append(
                    dict(
                        zip(
                            ["Group", "Sample ID", "Glu", "TG", "ALT"],
                            [group, len(rows) + 1, *values],
                        )
                    )
                )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Gtest.xlsx"
            pd.DataFrame(rows).to_excel(path, index=False)
            warnings = []
            result = compute_maha_model("reference", str(path), warnings)
            self.assertIsNotNone(result, warnings)
            self.assertEqual(len(result.individual_rows), 32)
            self.assertTrue(all(np.isfinite(row["NMD"]) for row in result.individual_rows))
            self.assertIn("HTD1801", result.htd_vs_pm)
            self.assertTrue(np.isfinite(result.htd_vs_pm["HTD1801"].g))
            write_part4_csv(
                str(Path(directory) / "distances.csv"),
                str(Path(directory) / "effects.csv"),
                SimpleNamespace(models=[result], pools={}),
            )
            with (Path(directory) / "effects.csv").open(encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            exported = float(
                next(row for row in rows if row["Comparison"] == "Model vs Control")["P_value"]
            )
            expected = two_sided_p(result.model_vs_control.g, result.model_vs_control.se)
            self.assertGreater(expected, 0)
            self.assertLess(expected, 0.00005)
            self.assertEqual(exported, expected)
