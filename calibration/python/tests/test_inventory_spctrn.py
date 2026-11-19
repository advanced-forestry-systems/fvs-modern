"""Regression checks for FIAJSP/SPCTRN resolution in inventory keyfiles."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import inventory_keyfile as ik


class InventorySpeciesResolutionTest(unittest.TestCase):
    def setUp(self):
        ik.reset_species_resolution_counts()
        ik.reset_ie_habitat_default_counts()
        self.stand = pd.DataFrame([{
            "age": 40, "site_index": 60, "forest_type": 161,
            "aspect": 0, "slope": 10, "elevft": 1200,
        }])

    def _trees(self, species):
        return pd.DataFrame([{
            "species": sp, "tree_count": 1.0, "diameter": 8.0,
            "ht": 40.0, "crratio": 50,
        } for sp in species])

    def test_sn_stage_two_rescues_stems_without_relabeling(self):
        trees = self._trees([68, 831, 110, 99999])
        before = ik.make_inventory_keyfile(self.stand, trees, "S1", "sn")
        self.assertEqual([line[14:17].strip() for line in before.splitlines()
                          if line.startswith("   1") and len(line) >= 35], ["SP"])
        ik.reset_species_resolution_counts()
        after = ik.make_inventory_keyfile(
            self.stand, trees, "S1", "sn", use_spctrn=True
        )
        records = [line[14:17].strip() for line in after.splitlines()
                   if line.startswith("   1") and len(line) >= 35]
        self.assertEqual(records, ["JU", "WK", "SP"])
        self.assertIn("     80301", after)
        counts = ik.get_species_resolution_counts("sn")
        self.assertEqual(counts["spctrn_stems"], 2)
        self.assertEqual(counts["fiajsp_stems"], 1)
        self.assertEqual(counts["dropped_stems"], 1)
        self.assertEqual(counts["dropped_by_spcd"][99999], 1)

    def test_ie_habitat_default_is_counted(self):
        ik.make_inventory_keyfile(
            self.stand, self._trees([202]), "IE1", "ie"
        )
        counts = ik.get_ie_habitat_default_summary()
        self.assertEqual(counts["total_ie_stands"], 1)
        self.assertEqual(counts["defaulted_stands"], 1)

    def test_unknown_variant_fails_closed(self):
        with self.assertRaises(KeyError):
            ik.make_inventory_keyfile(
                self.stand, self._trees([202]), "X", "unknown"
            )

    def test_spctrn_source_count_and_duplicate_fia_precedence(self):
        self.assertEqual(len(ik._spctrn_rows()), 562)
        self.assertEqual(ik.resolve_fvs_code(125, "ls"), ("RN", "fiajsp"))


if __name__ == "__main__":
    unittest.main()
