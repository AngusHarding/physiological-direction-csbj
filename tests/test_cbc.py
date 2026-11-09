"""Unit tests for CBC screening and derivation arithmetic."""

import unittest
import numpy as np
import pandas as pd


class CbcRulesTest(unittest.TestCase):
    def test_platelet_topcode_is_excluded(self):
        values = pd.Series([1.0, 975.0, 999.0, 999.9])
        self.assertListEqual(values.between(1, 999).tolist(), [True, True, True, False])

    def test_absolute_count(self):
        self.assertAlmostEqual(8.0 * 25.0 / 100, 2.0)

    def test_red_cell_indices(self):
        rbc, mcv, haemoglobin = 5.0, 90.0, 15.0
        hct = rbc * mcv / 10
        mch = haemoglobin * 10 / rbc
        mchc = haemoglobin * 100 / hct
        self.assertEqual((hct, mch), (45.0, 30.0))
        self.assertAlmostEqual(mchc, 33.333333333333336)

    def test_zero_lymphocyte_denominator_becomes_missing(self):
        ratio = pd.Series([1.0]) / pd.Series([0.0])
        ratio = ratio.replace([np.inf, -np.inf], np.nan)
        self.assertTrue(ratio.isna().iloc[0])


if __name__ == "__main__":
    unittest.main()
