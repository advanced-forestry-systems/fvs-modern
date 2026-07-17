#!/usr/bin/env python3
"""NE demonstration of the conus_greg consumer. Writes /tmp/greg_py_grid.csv."""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.config_loader import FvsConfigLoader
from sf_integration_dev.conus_greg_projector import GregProjector

SPECIES = {12: "Abies balsamea", 97: "Picea rubens", 316: "Acer rubrum",
           318: "Acer saccharum", 833: "Quercus rubra"}
DBHS = [6.0, 12.0, 20.0]
CR, BAL, CCFL, ELEV, CCH = 0.5, 90.0, 90.0, 1200.0, 0.4
EMT, TD = -28.8, 24.8

def ht_of(dbh):  # deterministic, shared with the R validator
    return 4.5 + 4.0 * dbh

L = FvsConfigLoader("ne", version="conus_greg")
gp = GregProjector(L, emt=EMT, td=TD)
rows = []
print("%-18s %4s %5s %8s %8s %8s %9s" %
      ("species", "SPCD", "DBH", "dg_ann", "hg_ann", "surv_ann", "dbh_5yr"))
for spcd, name in SPECIES.items():
    for dbh in DBHS:
        ht = ht_of(dbh)
        dg = gp.dg_annual(spcd, dbh, CR, ht, BAL, ELEV, emt=EMT)
        hg = gp.hg_annual(spcd, ht, CR, CCFL, CCH, ELEV, td=TD, emt=EMT)
        sv = gp.survival_annual(spcd, CR, CCH)
        d5 = gp.dg_increment(spcd, dbh, CR, ht, BAL, ELEV, emt=EMT, years=5)
        rows.append([spcd, name, dbh, dg, hg, sv, d5])
        print("%-18s %4d %5.0f %8.5f %8.5f %8.5f %9.5f" %
              (name, spcd, dbh, dg, hg, sv, d5))
with open("/tmp/greg_py_grid.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["SPCD", "name", "DBH", "dg_ann", "hg_ann", "surv_ann", "dbh_inc5"])
    w.writerows(rows)
print("\nwrote /tmp/greg_py_grid.csv (%d rows)" % len(rows))
