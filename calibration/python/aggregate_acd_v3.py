#!/usr/bin/env python3
"""
PERSEUS FVS-ACD v3 Aggregator
==============================
Aggregates per-batch FVS-ACD plot-level CSVs into state-total MMT trajectories
matching the methodology of perseus_aggregate.py / reaggregate_v3.py.

For each of the 6 PERSEUS cells, reads perseus_<cell>_batch*.csv files, dedups
by (PLOT, PROJ_YEAR, VARIANT, CONFIG), expands by approximate FIA EXPNS, and
writes perseus_<cell>_mmt_v3_acd.csv with state-total MMT_AGB / MMT (carbon) /
MEAN_AGB_TONS_AC / N_PLOTS columns.

Usage on Cardinal:
  source envs/libcbm/bin/activate || true   # only needs pandas
  python aggregate_acd_v3.py
"""
import glob
import os
import sys
import pandas as pd

ROOT = os.environ.get("FVS_PROJECT_ROOT",
                      "/users/PUOM0008/crsfaaron/fvs-modern")
OUT_DIR = os.path.join(ROOT, "calibration", "output", "perseus")

cells = ["noharvest_noclimate", "noharvest_rcp45", "noharvest_rcp85",
         "harvest_noclimate", "harvest_rcp45", "harvest_rcp85"]

MAINE_FOREST_ACRES = 17_600_000
CARBON_FRACTION = 0.5
TONS_TO_TONNES = 0.907185

for cell in cells:
    cell_dir = os.path.join(OUT_DIR, cell)
    batch_files = sorted(glob.glob(
        os.path.join(cell_dir, f"perseus_{cell}_batch*.csv")
    ))
    if not batch_files:
        print(f"  {cell}: no batch files found, skipping")
        continue

    dfs = []
    for f in batch_files:
        try:
            df = pd.read_csv(f)
            if "VARIANT" in df.columns:
                df = df[df["VARIANT"] == "ACD"]
                if len(df) > 0:
                    dfs.append(df)
        except Exception as e:
            print(f"    skip {f}: {e}")
    if not dfs:
        print(f"  {cell}: no ACD records found")
        continue

    combined = pd.concat(dfs, ignore_index=True)
    combined.drop_duplicates(
        subset=["PLOT", "PROJ_YEAR", "VARIANT", "CONFIG"], inplace=True
    )
    n_plots = combined["PLOT"].nunique()
    expns = MAINE_FOREST_ACRES / n_plots
    combined["EXPNS"] = expns

    # Aggregate to state-total per (VARIANT, CONFIG, PROJ_YEAR)
    g = combined.groupby(["VARIANT", "CONFIG", "PROJ_YEAR"], as_index=False)
    agg = g.apply(
        lambda d: pd.Series({
            "MMT_AGB": (d["AGB_TONS_AC"] * d["EXPNS"]).sum() * TONS_TO_TONNES / 1e6,
            "MMT":     (d["AGB_TONS_AC"] * d["EXPNS"]).sum() * TONS_TO_TONNES / 1e6 * CARBON_FRACTION,
            "MEAN_AGB_TONS_AC": d["AGB_TONS_AC"].mean(),
            "N_PLOTS": d["PLOT"].nunique(),
        }),
        include_groups=False,
    ).reset_index()
    agg["SCENARIO"] = cell

    out_path = os.path.join(cell_dir, f"perseus_{cell}_mmt_v3_acd.csv")
    agg.to_csv(out_path, index=False)

    # Brief summary
    cal = agg[(agg["VARIANT"] == "ACD") & (agg["CONFIG"] == "calibrated")].sort_values("PROJ_YEAR")
    if not cal.empty:
        v0 = cal.iloc[0]["MMT_AGB"]
        v50 = cal.loc[cal.PROJ_YEAR == 50, "MMT_AGB"]
        v100 = cal.loc[cal.PROJ_YEAR == 100, "MMT_AGB"]
        print(f"  {cell}: n_plots={n_plots}, EXPNS={expns:.0f}, "
              f"ACD_cal MMT_AGB y0={v0:.0f}, y50={v50.iloc[0] if len(v50)>0 else 'n/a':.0f}, "
              f"y100={v100.iloc[0] if len(v100)>0 else 'n/a':.0f}")

print("\nDone. v3_acd files written next to v3 NE files.")
