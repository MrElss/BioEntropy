import copy
import json
import unittest
import numpy as np
import pandas as pd
from localization import translate, presentation, localize_frame, sheet_title
from comparisons import compare_policies
from test_integrated_policy import fixture


class LocalizationTests(unittest.TestCase):
    def test_group_names_and_numeric_values_are_preserved(self):
        data = pd.DataFrame(
            {
                "Group": ["能量距离", "安慰剂"],
                "NMD": [1.25, 2.5],
                "Entropy_Shannon": [10.25, 11.5],
                "NMD_status": ["完整", "NA：全部指标缺失"],
            }
        )
        original = data.copy(deep=True)
        with presentation("en"):
            translated = localize_frame(data)
        pd.testing.assert_frame_equal(data, original)
        pd.testing.assert_series_equal(translated.Group, original.Group)
        np.testing.assert_array_equal(translated.NMD, data.NMD)
        np.testing.assert_array_equal(translated.Entropy_Shannon, data.Entropy_Shannon)
        self.assertEqual(translated.NMD_status.tolist(), ["Complete", "NA: all markers missing"])

    def test_raw_configuration_identifiers_are_preserved(self):
        metadata = {
            "features": ["指标1（单位）"],
            "groups": ["安慰剂"],
            "files": ["C:/数据/实验.xlsx"],
            "rows": [{"note": "6组"}],
        }
        data = pd.DataFrame(
            {
                "Setting": ["user_configuration", "nmd_features"],
                "Value": [json.dumps(metadata, ensure_ascii=False), '["指标1（单位）"]'],
            }
        )
        with presentation("en"):
            translated = localize_frame(data)
        result = json.loads(translated.Value.iloc[0])
        for key in ("features", "groups", "files"):
            self.assertEqual(result[key], metadata[key])
        self.assertEqual(result["rows"][0]["note"], "6 groups")
        self.assertEqual(translated.Value.iloc[1], data.Value.iloc[1])

    def test_localized_statistical_results(self):
        value = fixture()
        with presentation("zh"):
            expected = compare_policies(copy.deepcopy(value), 999, lambda text: None)
        with presentation("en"):
            actual = compare_policies(copy.deepcopy(value), 999, lambda text: None)
        self.assertEqual(expected, actual)

    def test_language_context_restoration(self):
        self.assertEqual(translate("显著性"), "显著性")
        with presentation("en"):
            self.assertEqual(sheet_title("图10_组值"), "Fig10_Groups")
            self.assertEqual(translate("显著性"), "Significance")
        self.assertEqual(translate("显著性"), "显著性")

    def test_source_label_preservation(self):
        with presentation("en", ["Shannon熵", "安慰剂"]):
            self.assertEqual(translate("Shannon熵"), "Shannon熵")
            self.assertEqual(translate("安慰剂"), "安慰剂")


if __name__ == "__main__":
    unittest.main()
