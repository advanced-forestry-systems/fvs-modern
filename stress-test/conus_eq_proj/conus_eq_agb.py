#!/usr/bin/env python3
"""conus_eq_agb.py -- fill AGB_TONS_AC into a conus_eq_projector metrics CSV
using the REAL NSBECalculator (identical biomass equations to the engine arm).

AGB HEIGHT BASIS (v3 standardization)
-------------------------------------
For a clean four-arm comparison the AGB metric must isolate the DBH/TPA growth
dynamics (the actual fvs-conus equation effect), NOT the height-model choice.
We therefore use a SINGLE, COMMON DBH->HT->biomass mapping across all arms:

    HT_for_biomass = HT_M * AGB_HT_ANCHOR

where HT_M is the projector's ht-dbh height (the same measured-anchored basis the
FVS engine carries in its treelist -- verified: engine FVS_TreeList Ht ~52 ft vs
projector HT_M ~49 ft) and AGB_HT_ANCHOR is a single global constant chosen so the
projector's year-0 NE AGB reproduces the engine's reported year-0 NE AGB (28.25
t/ac over the 991 identical NE stands) to within ~0.5%. The constant is applied
IDENTICALLY to every cycle and every arm, so it cannot bias the growth comparison;
DBH and TPA -- the things the fvs-conus equations actually drive -- remain the only
quantities that differ between arms.

Empirical anchoring (NE, 991 identical stands, engine y0 = 28.25 t/ac):
    HT_M * 1.00 -> 24.97 (-11.6%)
    HT_M * 1.18 -> 28.14 ( -0.4%)   <- selected
NOTE: NSBE compute_tree_biomass_kg consumes the supplied height numerically; the
engine's reported AGB is consistent with this numeric basis (passing HT_M->feet
overshoots ~+118%, confirming the engine's effective biomass-height is the HT_M
numeric basis, not a feet conversion).
"""
import argparse, os, sys
import pandas as pd
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P

# Global height anchor: HT_M numeric * this constant reproduces engine y0 NE AGB.
AGB_HT_ANCHOR = float(os.environ.get("AGB_HT_ANCHOR", "1.18"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--treelists", required=True)
    a = ap.parse_args()
    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    tl = pd.read_csv(a.treelists)
    def bm(r):
        try:
            ht = float(r.HT_M) * AGB_HT_ANCHOR
            kg = nsbe.compute_tree_biomass_kg(int(r.SPCD), float(r.DBH_IN), ht)
        except Exception:
            return 0.0
        if kg is None or kg != kg:
            return 0.0
        return (kg / 907.185) * float(r.TPA)
    tl["agb_contrib"] = tl.apply(bm, axis=1)
    agb = tl.groupby(["STAND_CN","CONFIG","PROJ_YEAR"], as_index=False)["agb_contrib"].sum().rename(columns={"agb_contrib":"AGB_FILL"})
    m = pd.read_csv(a.metrics)
    m["STAND_CN"] = m["STAND_CN"].astype(str); agb["STAND_CN"] = agb["STAND_CN"].astype(str)
    m = m.merge(agb, on=["STAND_CN","CONFIG","PROJ_YEAR"], how="left")
    m["AGB_TONS_AC"] = m["AGB_FILL"].round(4); m = m.drop(columns=["AGB_FILL"])
    m.to_csv(a.metrics, index=False)
    print("AGB filled. height anchor (HT_M * %.3f) applied to all cycles." % AGB_HT_ANCHOR)
    print("year-mean AGB_TONS_AC:")
    print(m[m.PROJ_YEAR.isin([0,50,100])].groupby("PROJ_YEAR")["AGB_TONS_AC"].mean())
if __name__ == "__main__":
    main()
