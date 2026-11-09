"""Tests for the independently reconstructed smoking domain."""

from pathlib import Path
import sys
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
import importlib.util  # public_code: script names start with a number, so load by path
_SPEC = importlib.util.spec_from_file_location("harmonise_smoking", ROOT / "4_lifestyle/02_harmonise_smoking.py")  # public_code
_MODULE = importlib.util.module_from_spec(_SPEC)  # public_code
_SPEC.loader.exec_module(_MODULE)  # public_code
classify = _MODULE.classify  # public_code


class SmokingTests(unittest.TestCase):
    def test_decision_table(self):
        example = pd.DataFrame({
            "SMQ020": [2, 1, 1, 1, 7, pd.NA],
            "SMQ040": [pd.NA, 1, 2, 3, pd.NA, pd.NA],
        })
        self.assertEqual(classify(example).tolist(), ["Never", "Current", "Current", "Former", pd.NA, pd.NA])

    def test_materialised_output(self):
        output = pd.read_csv(ROOT / "outputs/02_harmonised/smoking_status.csv")
        self.assertTrue(output.SEQN.is_unique)
        self.assertFalse(output.SEQN.isna().any())
        self.assertEqual(set(output.Smoking_status), {"Never", "Former", "Current"})


if __name__ == "__main__":
    unittest.main()

