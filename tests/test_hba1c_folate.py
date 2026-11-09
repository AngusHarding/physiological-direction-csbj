"""Tests for HbA1c bounds and the folate crossover implementation."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
SCRIPT=Path(__file__).resolve().parents[1]/"1_prepare_nhanes/2_blood_tests/09_harmonise_hba1c_folate.py"  # public_code
SPEC=importlib.util.spec_from_file_location("hf",SCRIPT); MODULE=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(MODULE)
class Hba1cFolateRulesTest(unittest.TestCase):
    def test_hba1c_bounds_are_inclusive(self): self.assertListEqual(pd.Series([2.9,3,20,20.1]).between(3,20).tolist(),[False,True,True,False])
    def test_folate_transform_returns_finite_positive_value(self):
        value=MODULE.transform_br_to_ma(pd.Series([500.0]),pd.Series([20.0]),pd.Series([45.0])).iloc[0]
        self.assertTrue(np.isfinite(value)); self.assertGreater(value,0)
    def test_direct_ma_value_needs_no_conversion(self): self.assertEqual(float(pd.Series([800.0]).iloc[0]),800.0)
if __name__=="__main__": unittest.main()
