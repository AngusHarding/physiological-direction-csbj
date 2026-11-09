from pathlib import Path
import sys
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
import importlib.util  # public_code: script names start with a number, so load by path
_SPEC = importlib.util.spec_from_file_location("harmonise_alcohol", ROOT / "4_lifestyle/04_harmonise_alcohol.py")  # public_code
_MODULE = importlib.util.module_from_spec(_SPEC)  # public_code
_SPEC.loader.exec_module(_MODULE)  # public_code
derive = _MODULE.derive  # public_code


class AlcoholTests(unittest.TestCase):
    def test_old_and_new_instrument_examples(self):
        example = pd.DataFrame({
            "SEQN": [1, 2, 3], "cycle": ["old", "old", "new"],
            "instrument_era": ["count_unit", "count_unit", "categorical_frequency"],
            "ALQ101": [1, 2, pd.NA], "ALQ110": [1, pd.NA, pd.NA],
            "ALQ120Q": [2, pd.NA, pd.NA], "ALQ120U": [1, pd.NA, pd.NA],
            "ALQ111": [pd.NA, pd.NA, 1], "ALQ121": [pd.NA, pd.NA, 5],
            "ALQ130": [3, pd.NA, 2],
        })
        result = derive(example)
        self.assertEqual(result.drinks_per_week.tolist(), [6.0, 0.0, 2.0])

    def test_materialised_output(self):
        output = pd.read_csv(ROOT / "outputs/02_harmonised/alcohol_drinks_per_week.csv")
        self.assertTrue(output.SEQN.is_unique)
        self.assertTrue(output.drinks_per_week.ge(0).all())


if __name__ == "__main__":
    unittest.main()

