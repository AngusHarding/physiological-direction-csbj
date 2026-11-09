from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class StandardBiochemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pd.read_csv(ROOT / "outputs" / "02_harmonised" / "standard_biochemistry.csv")
        cls.changes = pd.read_csv(ROOT / "outputs" / "02_harmonised" / "standard_biochemistry_change_log.csv")

    def test_certified_row_count(self):
        self.assertEqual(len(self.data), 63570)

    def test_unique_key_and_ten_cycles(self):
        self.assertFalse(self.data.SEQN.duplicated().any())
        self.assertEqual(self.data.cycle.nunique(), 10)

    def test_complete_positive_analytes(self):
        markers = self.data.columns.difference(["SEQN", "cycle"])
        self.assertFalse(self.data[markers].isna().any().any())
        self.assertTrue((self.data[markers] > 0).all().all())

    def test_all_14_calibrations_recorded(self):
        self.assertEqual(self.changes.rule_id.astype(str).str.startswith("BIO-CAL").sum(), 14)


if __name__ == "__main__":
    unittest.main()

