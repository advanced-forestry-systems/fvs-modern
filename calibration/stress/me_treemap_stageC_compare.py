#!/usr/bin/env python3
"""ME TreeMap pilot, Stage C: paint Maine AGB and compare FIADB vs TreeMap.

Joins the donor-plot FVS trajectories (Stage B) to the donor pixel areas
(Stage A) and aggregates Maine total above-ground biomass / carbon two ways,
isolating the AREA-EXPANSION choice (exactly like ycx_fiadb_vs_treemap.R, but
with FVS as the growth engine):

  TreeMap (spatial)   : sum_d  density_d(year) x area_ha_d
  FIADB  (uniform)    : mean_d density_d(year) x total_area_ha
                        (every donor plot weighted by the average area)

Both cover the same total Maine forest area; they diverge where a plot's biomass
correlates with how much area TreeMap assigns it. Output a comparison CSV + a
trajectory figure, per config (default vs calibrated).

Usage:
  python me_treemap_stageC_compare.py --traj me_donor_trajectories.csv \
      --donors me_treemap_donors.csv --start 2025 --out-dir me_pilot_out
"""
from __future__ import annotations
import argparse, os
import numpy as np, pandas as pd

AC_PER_HA = 2.4710538
TONS_TO_TG = 9.07185e-7       # short tons -> Tg (1 short ton = 0.907185 Mg)
C_FRACTION = 0.47


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True,
                    help="comma-separated trajectory CSV(s) "
                         "(default/calibrated + gompit)")
    ap.add_argument("--donors", required=True)
    ap.add_argument("--start", type=int, default=2025)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    tr = pd.concat([pd.read_csv(p, dtype={"PLT_CN": str})
                    for p in a.traj.split(",") if os.path.exists(p)],
                   ignore_index=True)
    dn = pd.read_csv(a.donors, dtype={"PLT_CN": str})[["PLT_CN", "area_ha"]]
    # one plot may back several TM_IDs -> sum area per PLT_CN
    area = dn.groupby("PLT_CN")["area_ha"].sum().reset_index()
    total_area_ha = area["area_ha"].sum()
    tr = tr.merge(area, on="PLT_CN", how="inner")
    tr["year"] = a.start + tr["PROJ_YEAR"].astype(int)
    # acres each plot represents
    tr["acres"] = tr["area_ha"] * AC_PER_HA

    rows = []
    for cfg, g in tr.groupby("CONFIG"):
        for yr, gy in g.groupby("year"):
            # TreeMap: biomass-weighted by actual pixel area
            tons_tm = float((gy["AGB_TONS_AC"] * gy["acres"]).sum())
            # FIADB: mean density x total area (uniform area per plot)
            mean_dens = float(gy["AGB_TONS_AC"].mean())
            tons_fia = mean_dens * (total_area_ha * AC_PER_HA)
            rows.append({
                "config": cfg, "year": int(yr),
                "agb_treemap_Tg": round(tons_tm * TONS_TO_TG, 3),
                "agb_fiadb_Tg":   round(tons_fia * TONS_TO_TG, 3),
                "agc_treemap_Tg": round(tons_tm * TONS_TO_TG * C_FRACTION, 3),
                "agc_fiadb_Tg":   round(tons_fia * TONS_TO_TG * C_FRACTION, 3),
                "n_plots": int(gy["PLT_CN"].nunique())})
    out = pd.DataFrame(rows).sort_values(["config", "year"])
    out["tm_over_fia"] = (out["agb_treemap_Tg"] / out["agb_fiadb_Tg"]).round(3)
    csv = os.path.join(a.out_dir, "me_fiadb_vs_treemap.csv")
    out.to_csv(csv, index=False)
    print(f"total ME forest area: {total_area_ha:,.0f} ha "
          f"({total_area_ha*AC_PER_HA:,.0f} ac)")
    print(out.to_string(index=False))

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5))
        col = {"default": "#999999", "calibrated": "#0072B2", "gompit": "#D55E00"}
        for cfg, g in out.groupby("config"):
            ax.plot(g["year"], g["agb_treemap_Tg"], "-", color=col.get(cfg, "k"),
                    label=f"{cfg} TreeMap")
            ax.plot(g["year"], g["agb_fiadb_Tg"], "--", color=col.get(cfg, "k"),
                    label=f"{cfg} FIADB")
        ax.set_xlabel("year"); ax.set_ylabel("Maine AGB (Tg dry)")
        ax.set_title("Maine 100-yr AGB: FVS x (TreeMap spatial vs FIADB uniform)")
        ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(os.path.join(a.out_dir, "me_fiadb_vs_treemap.png"), dpi=200)
        print("wrote me_fiadb_vs_treemap.png")
    except Exception as e:
        print("figure skipped:", e)


if __name__ == "__main__":
    main()
