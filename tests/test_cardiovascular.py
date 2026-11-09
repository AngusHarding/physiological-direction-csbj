"""Small synthetic tests for the decisions most likely to be misunderstood."""

import unittest
import numpy as np
import pandas as pd


class CardiovascularRulesTest(unittest.TestCase):
    def test_zero_is_removed_before_mean(self):
        readings = pd.DataFrame({"r2": [120.0], "r3": [0.0], "r4": [130.0]})
        readings = readings.mask(readings.eq(0))
        self.assertEqual(readings.mean(axis=1).iloc[0], 125.0)

    def test_available_readings_are_averaged(self):
        readings = pd.DataFrame({"r2": [120.0], "r3": [126.0], "r4": [np.nan]})
        self.assertEqual(readings.mean(axis=1).iloc[0], 123.0)

    def test_pulse_regularity_missing_stays_missing(self):
        values = pd.Series([1.0, 2.0, np.nan]).map({1.0: 0.0, 2.0: 1.0})
        self.assertEqual(values.iloc[0], 0.0)
        self.assertEqual(values.iloc[1], 1.0)
        self.assertTrue(pd.isna(values.iloc[2]))

    def test_error_bounds_are_inclusive(self):
        values = pd.Series([59.0, 60.0, 290.0, 291.0])
        self.assertListEqual(values.between(60, 290).tolist(), [False, True, True, False])


if __name__ == "__main__":
    unittest.main()
