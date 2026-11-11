import unittest
import numpy as np
from entropy_methods import shannon_animal, shannon_clinical, renyi2
from sample_information import component_scores, add_scores
from comparison_statistics import assignments, entropy_labelings, holm, permutation_p, compare
from display_coordinates import coordinates
from plots import aligned_table, display_table


class ScoreTests(unittest.TestCase):
    def test_animal_components_average_to_existing_estimate(self):
        x = np.random.default_rng(12).normal(size=(7, 3))
        scores, _, count = component_scores(x, 0.8, "animal", power=10, reps=3)
        h, _ = shannon_animal(x, 0.8, power=10, reps=3)
        self.assertAlmostEqual(float(scores.mean()), h, places=12)
        self.assertTrue((count == 3072).all())

    def test_clinical_quadrature_precision(self):
        x = np.random.default_rng(1).normal(size=(9, 3))
        scores, _, _ = component_scores(x, 0.8, "clinical")
        h, _ = shannon_clinical(x, 0.8)
        self.assertLess(abs(scores.mean() - h), 0.001)

    def test_shannon_inputs_and_missing_records(self):
        r = dict(
            name="test",
            mode="animal",
            _z=np.array([[0.0, 0.0], [1.0, 1.0], [np.nan, 2.0]]),
            summary=[dict(Group="g", Week=0, N_entropy=2, Bandwidth=0.8, Entropy_Shannon=2.0)],
            observations=[
                dict(Group="g", NMD=v, Included_entropy=i < 2) for i, v in enumerate([10, 20, 30])
            ],
        )
        original = r["summary"][0]["Entropy_Shannon"]
        add_scores(r, lambda s: None)
        self.assertTrue(np.isnan(r["observations"][2]["Shannon_information_score"]))
        self.assertEqual(r["summary"][0]["Entropy_Shannon"], original)
        first = r["observations"][0]
        self.assertAlmostEqual(
            first["Shannon_weighted_contribution"] * 2, first["Shannon_information_score"]
        )

    def test_affine_coordinates(self):
        r = dict(name="test", mode="clinical", summary=[], observations=[])
        for g, base in [("a", 10.0), ("b", 20.0)]:
            for w, change in [(0, 0.0), (6, -1.0)]:
                r["summary"].append(
                    dict(
                        Group=g,
                        Week=w,
                        NMD=base + change,
                        Entropy_Shannon=base + change,
                        Entropy_Renyi2=base + change,
                    )
                )
                r["observations"].append(
                    dict(Group=g, Week=w, NMD=base + 5.0, Shannon_information_score=base + 7.0)
                )
        aligned, _ = aligned_table([r])
        displays = display_table([r], aligned)
        groups, points = coordinates(r, "Entropy_Shannon", displays)
        self.assertEqual(groups[groups.Week.eq(0)].X.nunique(), 1)
        self.assertEqual(groups[groups.Week.eq(0)].X.iloc[0], 100)
        self.assertTrue((points.X > 100).all())
        self.assertTrue((points.Y > 100).all())
        _, renyi = coordinates(r, "Entropy_Renyi2", displays)
        self.assertFalse(renyi.Joint_plot_included.any())


class InferenceTests(unittest.TestCase):
    def test_random_labels_keep_trial_counts_and_patient_block(self):
        g = np.array([True, True, False, False, True, False])
        trial = np.array(["105"] * 4 + ["106"] * 2)
        a = assignments(g, trial, 999)
        np.testing.assert_array_equal(a[:, 0], g)
        self.assertTrue((a[:4].sum(0) == 2).all())
        self.assertTrue((a[4:].sum(0) == 1).all())
        same_patient_two_visits = np.array([0, 1, 2, 0, 1, 2])
        np.testing.assert_array_equal(
            a[same_patient_two_visits[:3]], a[same_patient_two_visits[3:]]
        )

    def test_entropy_permutation_is_correct_for_identical_kernels(self):
        x = np.zeros((6, 2))
        a = assignments([True] * 3 + [False] * 3, ["all"] * 6, 19)
        out = entropy_labelings(x, 0.7, a, 5)
        np.testing.assert_allclose(
            out["Entropy_Shannon"], np.log(2 * np.pi * np.e * 0.7**2), rtol=0, atol=1e-13
        )
        np.testing.assert_allclose(
            out["Entropy_Renyi2"], np.log(4 * np.pi * 0.7**2), rtol=0, atol=1e-13
        )

    def test_entropy_permutation_preserves_shift_and_scale(self):
        x = np.random.default_rng(9).normal(size=(8, 2))
        a = assignments([True] * 4 + [False] * 4, ["all"] * 8, 9)
        first = entropy_labelings(x, 0.8, a, 8)
        second = entropy_labelings(3 * x + 10, 2.4, a, 8)
        for key in first:
            np.testing.assert_allclose(second[key] - first[key], 2 * np.log(3), atol=1e-11)
        for gi in (0, 1):
            mask = a[:, 0] == 1 - gi
            self.assertAlmostEqual(
                first["Entropy_Renyi2"][0, gi], renyi2(x[mask], 0.8)[0], places=12
            )

    def test_plus_one_and_holm(self):
        self.assertEqual(permutation_p([10, 1, 2, 3])[0], 0.25)
        self.assertEqual(permutation_p([0, 0, 0, 0])[0], 1.0)
        np.testing.assert_allclose(holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
        self.assertEqual(holm([0.01, np.nan])[0], 0.02)

    def test_invalid_permutation_status(self):
        a = np.array([[1, 1], [1, 1], [0, 1], [0, 0]], float)
        values = entropy_labelings(np.zeros((4, 2)), 1.0, a, 5)
        self.assertTrue(np.isnan(values["Entropy_Shannon"]).all())
        self.assertTrue(np.isnan(permutation_p([1.0, np.nan, 2.0])[0]))

    def test_reference_and_disabled_inference(self):
        x = np.arange(12, dtype=float).reshape(6, 2) / 10
        r = dict(
            name="fixture",
            mode="animal",
            _z=x,
            roles=dict(treatment="A", comparator="B", reference="A"),
            observations=[
                dict(
                    Group="A" if i < 3 else "B",
                    Week=0,
                    Sample_ID=i,
                    NMD=float(i),
                    Included_entropy=True,
                )
                for i in range(6)
            ],
            summary=[
                dict(
                    Group=g,
                    Week=0,
                    NMD=1.0,
                    Entropy_Shannon=2.0,
                    Entropy_Renyi2=1.5,
                    N_NMD=3,
                    N_entropy=3,
                    Bandwidth=0.7,
                )
                for g in ["A", "B"]
            ],
        )
        rows = compare(r, 0, lambda s: None)
        self.assertTrue(all((np.isnan(row["P_raw"]) for row in rows)))
        rows = compare(r, 999, lambda s: None)
        nmd = next((row for row in rows if row["Metric"] == "NMD"))
        self.assertTrue(np.isnan(nmd["P_raw"]))
        self.assertIn("正常参照", nmd["Status"])


if __name__ == "__main__":
    unittest.main()
