from pathlib import Path
import sys
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
import importlib.util  # public_code: script names start with a number, so load by path
_SPEC = importlib.util.spec_from_file_location("harmonise_physical_activity", ROOT / "4_lifestyle/06_harmonise_physical_activity.py")  # public_code
_MODULE = importlib.util.module_from_spec(_SPEC)  # public_code
_SPEC.loader.exec_module(_MODULE)  # public_code
derive_early = _MODULE.derive_early  # public_code
derive_gpaq = _MODULE.derive_gpaq  # public_code


class PhysicalActivityTests(unittest.TestCase):
    def test_early_mapping(self):
        frame = pd.DataFrame({"SEQN": range(6), "instrument_era": "PAQ180",
                              "PAQ180": [1, 2, 3, 4, 7, pd.NA]})
        result = derive_early(frame).physically_active.tolist()
        self.assertEqual(result[:4], [0.0, 1.0, 1.0, 1.0])
        self.assertTrue(pd.isna(result[4]) and pd.isna(result[5]))

    def test_materialised_output(self):
        output = pd.read_csv(ROOT / "outputs/02_harmonised/physical_activity.csv")
        self.assertTrue(output.SEQN.is_unique)
        self.assertTrue(output.physically_active.dropna().isin([0, 1]).all())
        self.assertTrue(output.healthy_exercise.dropna().isin([0, 1]).all())

    def test_corrected_output_preserves_high_activity_respondents(self):
        historical = pd.read_csv(ROOT / "outputs/02_harmonised/physical_activity.csv")
        corrected = pd.read_csv(ROOT / "outputs/02_harmonised/physical_activity_inclusion_preserving.csv")
        self.assertEqual(len(corrected) - len(historical), 4260)
        self.assertLessEqual(corrected.met_minutes_week.max(), 10_000)
        self.assertTrue(corrected.SEQN.is_unique)
        flagged = corrected.met_extreme_above_10000_flag.eq(1)
        self.assertEqual(int(flagged.sum()), 4260)
        self.assertTrue(corrected.loc[flagged, "met_minutes_week_raw"].gt(10_000).all())


if __name__ == "__main__":
    unittest.main()
