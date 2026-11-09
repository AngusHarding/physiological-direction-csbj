"""Tests for anthropometric equations and partial-row retention."""

import unittest
import numpy as np
import pandas as pd


class AnthropometryRulesTest(unittest.TestCase):
    def test_bmi_uses_height_in_metres(self):
        self.assertAlmostEqual(70 / (175 / 100) ** 2, 22.857142857142858)

    def test_whtr_is_waist_over_height(self):
        self.assertEqual(85 / 170, 0.5)

    def test_partial_rows_are_retained(self):
        frame = pd.DataFrame({"weight": [70, np.nan], "height": [np.nan, np.nan], "waist": [np.nan, np.nan]})
        retained = frame.dropna(subset=list(frame), how="all")
        self.assertEqual(len(retained), 1)

    def test_bmi_and_whtr_require_different_numerators(self):
        weight, height, waist = np.nan, 170.0, 85.0
        self.assertTrue(np.isnan(weight / (height / 100) ** 2))
        self.assertEqual(waist / height, 0.5)


if __name__ == "__main__":
    unittest.main()
