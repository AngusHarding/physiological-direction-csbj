from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage


class SurveyDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pd.read_csv(ROOT / "outputs/01_imported/survey_design_1999_2018.csv")

    def test_unique_participants_and_ten_cycles(self):
        self.assertFalse(self.data.SEQN.duplicated().any())
        self.assertEqual(self.data.cycle.nunique(), 10)

    def test_pooled_weights_are_one_tenth(self):
        difference = (self.data.WTMEC20YR - self.data.WTMEC2YR / 10).abs().max()
        self.assertLess(difference, 1e-9)

    def test_design_fields_complete_for_positive_weight_records(self):
        positive = self.data.WTMEC20YR.gt(0)
        self.assertFalse(self.data.loc[positive, ["SDMVSTRA", "SDMVPSU", "survey_cluster"]].isna().any().any())


if __name__ == "__main__":
    unittest.main()
