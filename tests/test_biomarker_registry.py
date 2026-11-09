"""Guard the identity and roles of the manuscript's final 38 biomarkers."""

from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "1_prepare_nhanes"  # public_code: the preparation stage

# Keeping the expected identities explicit makes an accidental addition,
# deletion or role change fail loudly during routine repository testing.
EXPECTED_FINAL_38 = {
    "LBXSCA", "LBXSPH", "systolic", "diastolic", "BPXPLS",
    "LBDLYMNO", "LBDNENO", "LBDMONO", "LBDEONO", "LBDBANO", "LBXSGB",
    "LBXSBU", "LBXSUA", "LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSTB",
    "LBXSLDSI", "LBXSGTSI", "LBXSAPSI", "LBXGH", "WHtR", "LBXTC",
    "LBDHDD", "LBXSTR", "LBXPLTSI", "LBXMPSI", "LBXHGB", "LBXMCVSI",
    "LBXMC", "LBXRDW", "RBC_folate", "LBXSNASI", "LBXSKSI", "LBXSCLSI",
    "LBXSC3SI", "ACR", "eGFR",
}


class BiomarkerRegistryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = pd.read_csv(ROOT / "config/biomarker_registry.csv")
        cls.final = cls.registry.loc[cls.registry.final_38_member.eq("yes")]

    def test_registry_marker_names_are_unique(self):
        duplicates = self.registry.loc[self.registry.marker.duplicated(False), "marker"].tolist()
        self.assertEqual(duplicates, [])

    def test_final_panel_is_exactly_the_frozen_38(self):
        self.assertSetEqual(set(self.final.marker), EXPECTED_FINAL_38)
        self.assertEqual(len(self.final), 38)

    def test_every_final_marker_is_scored(self):
        self.assertTrue(self.final.complete_case_role.eq("scored").all())

    def test_derived_differentials_have_derivation_provenance(self):
        markers = {"LBDLYMNO", "LBDNENO", "LBDMONO", "LBDEONO", "LBDBANO"}
        rows = self.registry[self.registry.marker.isin(markers)]
        self.assertEqual(len(rows), 5)
        self.assertTrue(rows.measured_or_derived.eq("derived").all())
        self.assertTrue(rows.source_component.eq("cbc_derived").all())


if __name__ == "__main__":
    unittest.main()
