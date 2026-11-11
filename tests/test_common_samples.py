import unittest

import numpy as np
import pandas as pd

from common_samples import common_records


class SharedSamplesTest(unittest.TestCase):
    def test_animal_partial_panels_and_singleton_group(self):
        dist, keep = common_records(
            [1, 2, 3, 4, 5],
            pd.Series([True, False, True, False, True]),
            ["A", "A", "A", "B", "B"],
            0,
        )
        np.testing.assert_equal(keep, [True, False, True, False, False])
        np.testing.assert_allclose(dist, [1, np.nan, 3, np.nan, np.nan])

    def test_missing_distance_is_local_to_visit(self):
        dist, keep = common_records(
            [[4, np.nan], [2, 3], [1, 2]],
            np.ones((3, 2), bool),
            np.array(["A", "A", "A"])[:, None],
            [0, 6],
        )
        np.testing.assert_equal(keep, [[True, False], [True, True], [True, True]])
        np.testing.assert_allclose(dist, [[4, np.nan], [2, 3], [1, 2]])
