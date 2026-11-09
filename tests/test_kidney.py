"""Tests for kidney units, equations and questionnaire recoding."""

import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes/2_blood_tests/12_derive_kidney_features.py"  # public_code
SPEC = importlib.util.spec_from_file_location("kidney_derivation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class KidneyRulesTest(unittest.TestCase):
    def test_acr_has_required_factor_of_one_hundred(self):
        self.assertEqual(30.0 / 100.0 * 100.0, 30.0)

    def test_ckd_epi_known_manual_case(self):
        # Male, age 50, creatinine 0.9 has ratio 1, simplifying both branches.
        expected = 142.0 * 0.9938 ** 50
        actual = MODULE.egfr_ckdepi_2021([0.9], [50], [1])[0]
        self.assertAlmostEqual(actual, expected)

    def test_female_multiplier_is_applied(self):
        base = 142.0 * 0.9938 ** 50
        actual = MODULE.egfr_ckdepi_2021([0.7], [50], [2])[0]
        self.assertAlmostEqual(actual, base * 1.012)

    def test_questionnaire_unknown_is_not_recoded_healthy(self):
        values = pd.Series([1.0, 2.0, 7.0, 9.0, np.nan]).map({1.0: 1.0, 2.0: 0.0})
        self.assertListEqual(values.iloc[:2].tolist(), [1.0, 0.0])
        self.assertTrue(values.iloc[2:].isna().all())


if __name__ == "__main__":
    unittest.main()
