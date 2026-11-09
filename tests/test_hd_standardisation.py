"""Formula and structure locks for healthy-reference HD standardisation."""

from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage
Z_PATH = ROOT / "outputs/09_standardised/final_38_panel_transformed_standardised_hd.csv"
PARAM_PATH = ROOT / "outputs/09_standardised/hd_reference_standardisation_parameters.csv"
Y_PATH = ROOT / "outputs/06_transform/final_38_panel_transformed_unstandardised_v2.csv"


class HDStandardisationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z = pd.read_csv(Z_PATH)
        cls.y = pd.read_csv(Y_PATH)
        cls.params = pd.read_csv(PARAM_PATH)
        cls.markers = pd.read_csv(ROOT / "config/final_38_transformations_v2.csv").marker.tolist()

    def test_participant_and_marker_identity(self):
        self.assertEqual(self.z.shape, (29053, 42))
        self.assertListEqual(list(self.z.columns),
                             ["SEQN", "RIAGENDR", "RIDAGEYR", "is_reference"] + self.markers)
        pd.testing.assert_series_equal(self.z.SEQN, self.y.SEQN)
        self.assertFalse(self.z.SEQN.duplicated().any())

    def test_reference_counts_and_parameter_rows(self):
        self.assertEqual(int(self.z.is_reference.sum()), 1238)
        self.assertEqual(len(self.params), 76)
        self.assertTrue(self.params.groupby("sex").reference_n.first().eq(
            pd.Series({"female": 691, "male": 547})).all())

    def test_reference_coordinates_are_mean_zero_sd_one(self):
        for sex_code in (1.0, 2.0):
            ref = self.z.RIAGENDR.eq(sex_code) & self.z.is_reference.astype(bool)
            self.assertLess(float(self.z.loc[ref, self.markers].mean().abs().max()), 1e-12)
            self.assertLess(float((self.z.loc[ref, self.markers].std(ddof=1) - 1).abs().max()), 1e-12)

    def test_saved_values_reproduce_parameters(self):
        lookup = self.params.set_index(["marker", "sex_code"])
        for marker in self.markers:
            for sex_code in (1, 2):
                mask = self.y.RIAGENDR.eq(sex_code)
                row = lookup.loc[(marker, sex_code)]
                expected = ((self.y.loc[mask, marker] - row.reference_mean_transformed)
                            / row.reference_sd_transformed)
                np.testing.assert_allclose(self.z.loc[mask, marker], expected,
                                           rtol=1e-12, atol=1e-12)

    def test_no_missing_or_infinite_coordinates(self):
        self.assertTrue(np.isfinite(self.z[self.markers].to_numpy()).all())


if __name__ == "__main__":
    unittest.main()
