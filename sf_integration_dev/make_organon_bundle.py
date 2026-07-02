#!/usr/bin/env python3
"""Generate a structurally-valid ORGANON-form DG bundle (arm 1 scaffold).

The ORGANON diameter-growth form (Hann/ORGANON DG_NWO family):
  ln(dDBH) = a0 + a1*ln(DBH) + a2*DBH^2 + a3*ln_cr_adj
           + a4*ln(SI) + a5*(BAL/ln(DBH+5)) + a6*sqrt(SBA)
Per-species multiplicative calibration ratios (ORGANON CALIB) start at 1.0.
Coefficients are NOMINAL placeholders pending the CONUS-wide Bayesian ORGANON
recalibration and Greg's climate-only DG form. Tagged accordingly.
"""
import json, sys

# nominal ORGANON DG_NWO-family coefficients (documented placeholder ranges)
fixed = {"a0": -3.20, "a1": 0.42, "a2": -0.00035, "a3": 0.98,
         "a4": 0.55, "a5": -0.015, "a6": -0.045}

# a small representative species set with neutral calibration (CALIB=1.0);
# the recalibration fit populates per-species calib ratios.
spcds = [int(x) for x in sys.argv[1:]] or [202, 17, 122, 93, 15]
bundle = {
    "model": "organon_dg",
    "form": ("ln(dDBH) = a0 + a1*ln_dbh + a2*dbh_sq + a3*ln_cr_adj "
             "+ a4*ln_si + a5*bal_over_lndbh + a6*sqrt_sba"),
    "fixed_effects": {"param": list(fixed.keys()),
                      "mean": list(fixed.values())},
    "species": {"SPCD": spcds, "calib": [1.0] * len(spcds)},
    "status": "scaffold_pending_fit",
    "notes": ("ORGANON-form (arm 1). Coefficients are NOMINAL placeholders. "
              "Backfill via 62c once the CONUS-wide Bayesian ORGANON "
              "recalibration and Greg's climate-only DG form land. Selected "
              "with version=conus_organon; the shared categories_conus_mod "
              "modifier layer applies on top.")
}
json.dump(bundle, open("/users/PUOM0008/crsfaaron/dg_organon_bundle.json", "w"),
          indent=2)
print("wrote dg_organon_bundle.json with", len(spcds), "species")
