import copy
import unittest
import numpy as np
from entropy_methods import shannon_clinical
from comparisons import compare_policies, shared_evaluations
import comparison_statistics as main_statistics
import visit_statistics as raw
from analysis_outputs import annotation
from weekly_state import coordinates as raw_coordinates
from display_coordinates import coordinates as aligned_coordinates
from plots import aligned_table, display_table


def fixture():
    rng = np.random.default_rng(518)
    rows, points, vectors = ([], [], [])
    for week in (0, 6):
        for arm, group in enumerate(("A", "B")):
            x = rng.normal(size=(4, 2)) * (1.3 if arm else 0.9) + week / 20
            distances = np.arange(4.0) + (8 if arm else 10) - week / (4 if arm else 2)
            h = shannon_clinical(x, 0.8)[0]
            rows.append(
                dict(
                    Group=group,
                    Week=week,
                    NMD=float(distances.mean()),
                    Entropy_Shannon=h,
                    N_NMD=4,
                    N_entropy=4,
                    Bandwidth=0.8,
                )
            )
            for i in range(4):
                points.append(
                    dict(
                        Group=group,
                        Week=week,
                        Trial="T",
                        USUBJID=group + str(i),
                        ID_verified=True,
                        Included_entropy=True,
                        NMD=distances[i],
                    )
                )
            vectors.extend(x)
    return dict(
        name="fixture",
        mode="clinical",
        selected_method="Shannon",
        summary=rows,
        observations=points,
        _z=np.array(vectors),
        roles=dict(treatment="A", comparator="B", reference="external"),
    )


class IntegratedPolicyTests(unittest.TestCase):
    def test_statistical_policy_dispatch(self):
        base = fixture()
        expected_main = main_statistics.compare(copy.deepcopy(base), 999, lambda s: None)
        expected_visit = raw.compare(copy.deepcopy(base), 999, lambda s: None)
        result = copy.deepcopy(base)
        combined = compare_policies(result, 999, lambda s: None)
        for source, expected, actual in [
            ("常规", expected_main, result["comparisons"]),
            ("分周原值", expected_visit, result["comparisons_visit"]),
        ]:
            self.assertEqual(len(expected), len(actual))
            for a, b in zip(expected, actual):
                for key in [
                    "Metric",
                    "Week",
                    "Contrast",
                    "P_raw",
                    "P_Holm",
                    "Annotation",
                    "Caution",
                    "Nodes_per_kernel",
                ]:
                    self.assertEqual(a[key], b[key], (source, key))
                self.assertEqual(b["Inference_source"], source)
        self.assertEqual(len(combined), len(expected_main) + len(expected_visit))
        self.assertGreater(result["inference_settings"]["shared_kernel_evaluations"]["reused"], 0)
        self.assertEqual(base["summary"], result["summary"])

    def test_annotation_policy_selection(self):
        r = dict(
            mode="clinical",
            inference_settings={"permutations": 999},
            comparisons=[
                dict(Metric="Entropy_Shannon", Week=6, Contrast="raw", Annotation="*"),
                dict(Metric="Entropy_Shannon", Week=6, Contrast="change", Annotation="**"),
            ],
            comparisons_visit=[
                dict(Metric="Entropy_Shannon", Week=6, Contrast="raw", Annotation="ns")
            ],
        )
        self.assertEqual(annotation(r, "Entropy_Shannon", 6, False), "*")
        self.assertEqual(annotation(r, "Entropy_Shannon", 6, True), "**")
        self.assertEqual(annotation(r, "Entropy_Shannon", 6, False, source="分周原值"), "ns")
        self.assertEqual(annotation(r, "Entropy_Shannon", 0, True), "0W对齐")

    def test_kernel_cache_is_exact_and_distinguishes_bandwidth_and_labels(self):
        x = np.random.default_rng(4).normal(size=(6, 2))
        labels = main_statistics.assignments([True] * 3 + [False] * 3, ["T"] * 6, 9)
        expected = main_statistics.entropy_labelings(
            x, 0.8, labels, 5, selected=["Entropy_Shannon"]
        )
        with shared_evaluations() as stats:
            first = main_statistics.entropy_labelings(
                x, 0.8, labels, 5, selected=["Entropy_Shannon"]
            )
            second = raw.entropy_labelings(
                x.copy(), 0.8, labels.copy(), 5, selected=["Entropy_Shannon"]
            )
            raw.entropy_labelings(x, 0.9, labels, 5, selected=["Entropy_Shannon"])
            raw.entropy_labelings(x, 0.8, 1 - labels, 5, selected=["Entropy_Shannon"])
        np.testing.assert_array_equal(expected["Entropy_Shannon"], first["Entropy_Shannon"])
        np.testing.assert_array_equal(first["Entropy_Shannon"], second["Entropy_Shannon"])
        self.assertEqual(stats, dict(evaluations=3, reused=1))

    def test_clinical_joint_points_are_raw_while_trend_policy_stays_aligned(self):
        r = dict(
            name="x",
            mode="clinical",
            selected_method="Shannon",
            summary=[
                dict(Group="A", Week=0, NMD=10.0, Entropy_Shannon=2.0),
                dict(Group="A", Week=6, NMD=6.0, Entropy_Shannon=1.0),
                dict(Group="B", Week=0, NMD=20.0, Entropy_Shannon=5.0),
                dict(Group="B", Week=6, NMD=19.0, Entropy_Shannon=4.0),
            ],
            observations=[
                dict(Group="A", Week=0, NMD=300.0, Shannon_information_score=-40.0),
                dict(Group="B", Week=0, NMD=2.0, Shannon_information_score=400.0),
            ],
        )
        aligned, _ = aligned_table([r])
        display = display_table([r], aligned)
        groups, points = raw_coordinates(r, "Entropy_Shannon", display)
        np.testing.assert_array_equal(groups.X, [2.0, 1.0, 5.0, 4.0])
        np.testing.assert_array_equal(groups.Y, [10.0, 6.0, 20.0, 19.0])
        np.testing.assert_array_equal(points.X, [-40.0, 400.0])
        old, _ = aligned_coordinates(r, "Entropy_Shannon", display)
        self.assertEqual(old[old.Week.eq(0)].X.nunique(), 1)

    def test_animal_uses_main_comparisons(self):
        r = fixture()
        r["mode"] = "animal"
        r["summary"] = [x for x in r["summary"] if x["Week"] == 0]
        r["observations"] = [x for x in r["observations"] if x["Week"] == 0]
        r["_z"] = r["_z"][:8]
        compare_policies(r, 0, lambda s: None)
        self.assertEqual(r["comparisons_visit"], [])
        self.assertEqual({x["Inference_source"] for x in r["comparisons_all"]}, {"常规"})
        self.assertEqual({x["Contrast"] for x in r["comparisons_all"]}, {"raw"})


if __name__ == "__main__":
    unittest.main()
