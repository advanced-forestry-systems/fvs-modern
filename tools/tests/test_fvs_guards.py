"""Regression checks for diagnostic FVS input and run guards."""
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import fvs_input_guards as ig
import fvs_run_guards as rg


class InputGuardTest(unittest.TestCase):
    def test_ne_species_table_matches_current_config(self):
        spec = json.loads((ROOT / "config/ne.json").read_text())["categories"]["species_definitions"]
        actual = tuple((int(f.strip()) if f.strip() else None, a.strip())
                       for f, a in zip(spec["FIAJSP"], spec["JSP"]))
        self.assertEqual(actual, ig._NE_TABLE)

    def test_species_and_stand_ids(self):
        self.assertEqual(ig.classify_spcd(97)["cls"], ig.VALID)
        self.assertEqual(ig.classify_spcd(263)["cls"], ig.OUT_OF_REGION)
        self.assertEqual(ig.classify_spcd(9999)["cls"], ig.INVALID)
        with self.assertRaises(ig.SpeciesCodeError):
            ig.validate_spcd(9999)
        ids = ig.build_standid_map(["10000001", "20000001"], "ne")
        self.assertEqual(ids.n_distinct(), 2)
        self.assertEqual(ig.legacy_standid("10000001"), ig.legacy_standid("20000001"))


class RunGuardTest(unittest.TestCase):
    def test_cycle_limit_and_fraction(self):
        self.assertEqual(rg.validate_numcycle("40.0", maxcyc=40), 40)
        for bad in (41, 2.5, "10.5", "bad", True, 0):
            with self.subTest(bad=bad), self.assertRaises(rg.FVSKeywordError):
                rg.validate_numcycle(bad, maxcyc=40)

    def test_keyword_and_listing(self):
        with tempfile.TemporaryDirectory() as d:
            key = Path(d) / "stand.key"
            listing = Path(d) / "stand.out"
            key.write_text("TIMEINT 0 1\nNUMCYCLE 10.0\n")
            self.assertEqual(rg.validate_keyword_file(str(key), maxcyc=40)["numcycle"], 10)
            for bad in ("NUMCYCLE 10.5\n", "NUMCYCLE\n", "NUMCYCLE 41\n"):
                with self.subTest(bad=bad), self.assertRaises(rg.FVSKeywordError):
                    key.write_text(bad)
                    rg.validate_keyword_file(str(key), maxcyc=40)
            listing.write_text("Normal completion\n")
            self.assertTrue(rg.check_run(10, str(listing), expected_cycles=10, summary_rows=11)["ok"])
            with self.assertRaises(rg.FVSRunError):
                rg.check_run(20, str(listing), expected_cycles=10, summary_rows=2)
            listing.write_text("FVS04 ERROR: KEYWORD IGNORED\n")
            with self.assertRaises(rg.FVSRunError):
                rg.check_run(10, str(listing), expected_cycles=10, summary_rows=11)
            with self.assertRaises(rg.FVSRunError):
                rg.check_run(10, str(Path(d) / "missing.out"))

    def test_summary_read_does_not_create_missing_database(self):
        with tempfile.TemporaryDirectory() as d:
            missing = Path(d) / "missing.db"
            self.assertIsNone(rg.summary_row_count(str(missing)))
            self.assertFalse(missing.exists())
            db = Path(d) / "fvsout.db"
            with sqlite3.connect(db) as con:
                con.execute("CREATE TABLE FVS_Summary2 (Year INTEGER)")
                con.executemany("INSERT INTO FVS_Summary2 VALUES (?)", [(2000,), (2005,)])
            self.assertEqual(rg.summary_row_count(str(db)), 2)


if __name__ == "__main__":
    unittest.main()
