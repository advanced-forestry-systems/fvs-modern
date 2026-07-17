#!/usr/bin/env python3
"""fvs_treemap_conus_compare.py -- CONUS FVS x (TreeMap spatial vs FIADB uniform).

Scales the Maine TreeMap pilot to CONUS. The v2 FVS campaign ran every FIA plot
(keyed STAND_CN = PLT_CN). TreeMap2022 imputes each FIA plot across the CONUS
forest landscape; its VAT gives the pixel COUNT per PLT_CN, hence the actual area
each plot represents (plt_area_treemap.csv, no raster scan). We expand the FVS
carbon two ways, isolating the AREA choice (exactly the ME pilot, CONUS-wide):

  TreeMap (spatial) : sum_plot  density_plot(year) x area_ha_plot
  FIADB  (uniform)  : mean_plot density_plot(year) x total_area     (same total
                      area, distributed uniformly across plots)

at state scale and by forest-type stratum (the varying-scale cut), for years
2030/2075/2125, default reserve engine. Also cross-checks the FVS 2030 TreeMap
carbon against TreeMap's own imputed live carbon (CARBON_L) for the same area.

Usage:
  python3 fvs_treemap_conus_compare.py --campaign out_fvs_v2 \
     --areas plt_area_treemap.csv --config default --out treemap_conus
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd

TONS_AC_TO_MGHA = 2.241702         # AGB tons/ac -> Mg/ha
C_FRACTION = 0.47
TONS_C_AC_TO_MGHA = 2.241702       # CARBON_L tons C/ac -> Mg C/ha
YEARS = {5: 2030, 50: 2075, 100: 2125}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--areas", required=True)
    ap.add_argument("--config", default="default")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    area = pd.read_csv(a.areas, dtype={"PLT_CN": str})
    area["PLT_CN"] = area["PLT_CN"].str.replace(r"\.0$", "", regex=True)

    # load campaign reserve densities at the 3 years
    rows = []
    for f in glob.glob(os.path.join(a.campaign, "conus_*.csv")):
        d = pd.read_csv(f, usecols=["STAND_CN", "STATE", "CONFIG",
                                    "PROJ_YEAR", "AGB_TONS_AC"])
        d = d[(d.CONFIG == a.config) & (d.PROJ_YEAR.isin(YEARS))]
        if len(d):
            rows.append(d)
    camp = pd.concat(rows, ignore_index=True)
    camp["PLT_CN"] = camp["STAND_CN"].astype(str).str.replace(r"\.0$", "",
                                                              regex=True)
    camp["agc_MgC_ha"] = camp.AGB_TONS_AC * TONS_AC_TO_MGHA * C_FRACTION
    camp = camp.merge(area[["PLT_CN", "area_ha", "tm_fortyp", "tm_carbon_L"]],
                      on="PLT_CN", how="inner")
    print(f"{camp.PLT_CN.nunique()} plots joined to TreeMap area "
          f"({camp.area_ha.sum()/3/1e6:.0f} Mha across 3 yrs)")

    # ---- state scale ----
    def expand(g):
        out = {}
        for py, yr in YEARS.items():
            gy = g[g.PROJ_YEAR == py]
            if not len(gy):
                continue
            tot_area = gy.area_ha.sum()
            tm = float((gy.agc_MgC_ha * gy.area_ha).sum()) / 1e6      # Tg C
            fia = float(gy.agc_MgC_ha.mean()) * tot_area / 1e6        # Tg C
            out[yr] = (tm, fia)
        return out

    st_rows = []
    for st, g in camp.groupby("STATE"):
        e = expand(g)
        for yr, (tm, fia) in e.items():
            st_rows.append({"scale": "state", "key": st, "year": yr,
                            "treemap_TgC": round(tm, 3), "fiadb_TgC": round(fia, 3),
                            "tm_over_fia": round(tm / fia, 3) if fia else None})
    # ---- CONUS total + forest-type strata ----
    e = expand(camp)
    for yr, (tm, fia) in e.items():
        st_rows.append({"scale": "CONUS", "key": "CONUS", "year": yr,
                        "treemap_TgC": round(tm, 3), "fiadb_TgC": round(fia, 3),
                        "tm_over_fia": round(tm / fia, 3) if fia else None})
    for ft, g in camp.groupby("tm_fortyp"):
        e = expand(g)
        for yr, (tm, fia) in e.items():
            st_rows.append({"scale": "fortyp", "key": int(ft), "year": yr,
                            "treemap_TgC": round(tm, 3), "fiadb_TgC": round(fia, 3),
                            "tm_over_fia": round(tm / fia, 3) if fia else None})
    out = pd.DataFrame(st_rows)
    out.to_csv(os.path.join(a.out, "fvs_treemap_vs_fiadb.csv"), index=False)

    # ---- cross-check: FVS 2030 TreeMap carbon vs TreeMap's own CARBON_L ----
    g0 = camp[camp.PROJ_YEAR == 5].drop_duplicates("PLT_CN")
    fvs_tm = float((g0.agc_MgC_ha * g0.area_ha).sum()) / 1e6
    tmc = g0.tm_carbon_L * TONS_C_AC_TO_MGHA
    tm_own = float((tmc * g0.area_ha).sum()) / 1e6
    print(f"CONUS 2030 live carbon: FVS-on-TreeMap {fvs_tm:.0f} TgC  vs  "
          f"TreeMap-native {tm_own:.0f} TgC  (ratio {fvs_tm/tm_own:.2f})")

    conus = out[out.scale == "CONUS"]
    print("\nCONUS FVS reserve carbon (Tg C):")
    print(conus[["year", "treemap_TgC", "fiadb_TgC", "tm_over_fia"]]
          .to_string(index=False))
    print(f"\nwrote {a.out}/fvs_treemap_vs_fiadb.csv "
          f"({out.scale.eq('state').sum()} state rows, "
          f"{out.scale.eq('fortyp').sum()} forest-type rows)")


if __name__ == "__main__":
    main()
