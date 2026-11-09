import importlib.util
from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


def load_module():
    path = ROOT / "5_health_conditions/02_harmonise_self_reported_disease.py"  # public_code
    spec = importlib.util.spec_from_file_location("disease", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class DiseaseTests(unittest.TestCase):
    def test_binary_preserves_unknown_as_missing(self):
        module = load_module()
        result = module.binary(pd.Series([1, 2, 3, 7, 9, float("nan")]))
        self.assertEqual(result.iloc[0], 1)
        self.assertEqual(result.iloc[1], 0)
        self.assertTrue(result.iloc[2:].isna().all())

    def test_primary_set_has_13_diseases_and_excludes_transfusion(self):
        module = load_module()
        self.assertEqual(len(module.PRIMARY_13), 13)
        self.assertNotIn("blood_transfusion", module.PRIMARY_13)

    def test_output_does_not_impute_unasked_osteoporosis(self):
        flags = pd.read_csv(ROOT / "outputs/02_harmonised/self_reported_disease_flags.csv")
        unfielded = flags.cycle.isin(["2011-2012", "2015-2016"])
        self.assertFalse(flags.loc[unfielded, "osteoporosis_observed"].any())
        self.assertTrue(flags.loc[unfielded, "osteoporosis"].isna().all())


if __name__ == "__main__":
    unittest.main()
