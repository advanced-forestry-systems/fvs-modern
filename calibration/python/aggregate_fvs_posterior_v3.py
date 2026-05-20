#!/usr/bin/env python3
"""
PERSEUS FVS Posterior v3 Aggregator
====================================

Reads per-batch CSVs from factorial_v3_posterior/<cell>/<cell>/perseus_<cell>_batch*.csv,
extracts posterior draw rows (CONFIG=='draw'), computes state-total MMT_AGB per
draw per year, then aggregates to q10/q25/median/q75/q90 bands per cell per
variant.

State-total: MMT_AGB = mean(AGB_TONS_AC across surviving plots) * 17.6M acres
* 0.907 / 1e6 (per-year EXPNS, matches v3 + v3_acd methodology).

Usage on Cardinal:
  python3 aggregate_fvs_posterior_v3.py
"""
import glob
import os
import sys
import pandas as pd
import numpy as np

ROOT = os.environ.get("FVS_PROJECT_ROOT",
                      "/users/PUOM0008/crsfaaron/fvs-modern")
OUT_BASE = os.path.join(ROOT, "calibration", "output", "perseus", "factorial_v3_posterior")

cells = ["noharvest_noclimate", "noharvest_rcp45", "noharvest_rcp85",
         "harvest_noclimate", "harvest_rcp45", "harvest_rcp85"]

MAINE_FOREST_ACRES = 17_600_000
TONS_TO_TONNES = 0.907185
CARBON_FRACTION = 0.5
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
QLABELS = ["q10", "q25", "median", "q75", "q90"]


def aggregate_cell(cell):
    """Aggregate one cell to per-variant per-year posterior bands."""
    cell_dir = os.path.join(OUT_BASE, cell, cell)
    if not os.path.isdir(cell_dir):
        cell_dir = os.path.join(OUT_BASE, cell)
    files = sorted(glob.glob(os.path.join(cell_dir, f"perseus_{cell}_batch*.csv")))
    files = [f for f in files if "_mmt_" not in f]
    if not files:
        print(f"  {cell}: no per-batch CSVs found at {cell_dir}")
        return None
    print(f"  {cell}: reading {len(files)} batch files ...")
    dfs = []
    for f in files:
        try:
            d = pd.read_csv(f)
            dfs.append(d)
        except Exception as e:
            print(f"    skip {f}: {e}")
    df = pd.concat(dfs, ignore_index=True)

    # Posterior draw rows have CONFIG='draw' and a DRAW column
    draws = df[df["CONFIG"] == "draw"].copy()
    if draws.empty:
        print(f"    no posterior draw rows found in {cell}")
        return None
    if "DRAW" not in draws.columns:
        print(f"    no DRAW column found in {cell}")
        return None

    # Per-year EXPNS based on surviving plots at each (variant, year, draw)
    rows = []
    g = draws.groupby(["VARIANT", "PROJ_YEAR", "DRAW"], as_index=False)
    print(f"    aggregating {len(g)} groups ...", flush=True)
    for (var, py, draw), grp in g:
        n = grp["PLOT"].nunique()
        if n == 0:
            continue
        mean_agb = grp["AGB_TONS_AC"].mean()
        total_tons = mean_agb * MAINE_FOREST_ACRES
        mmt_agb = total_tons * TONS_TO_TONNES / 1e6
        rows.append({
            "VARIANT": var, "PROJ_YEAR": int(py), "DRAW": int(draw),
            "MMT_AGB": mmt_agb,
        })
    state = pd.DataFrame(rows)
    print(f"    {cell}: state-total draws computed: {len(state)} (variant x year x draw)")

    # Quantile bands per variant per year
    bands_rows = []
    for (var, py), grp in state.groupby(["VARIANT", "PROJ_YEAR"]):
        vals = grp["MMT_AGB"].values
        n_draws = len(vals)
        qs = np.quantile(vals, QUANTILES)
        for ql, qv in zip(QLABELS, qs):
            bands_rows.append({"VARIANT": var, "PROJ_YEAR": int(py),
                               "STAT": ql, "VALUE": float(qv),
                               "N_DRAWS": int(n_draws), "SCENARIO": cell})
        bands_rows.append({"VARIANT": var, "PROJ_YEAR": int(py),
                           "STAT": "mean", "VALUE": float(np.mean(vals)),
                           "N_DRAWS": int(n_draws), "SCENARIO": cell})
        bands_rows.append({"VARIANT": var, "PROJ_YEAR": int(py),
                           "STAT": "sd", "VALUE": float(np.std(vals, ddof=1)),
                           "N_DRAWS": int(n_draws), "SCENARIO": cell})
    bands = pd.DataFrame(bands_rows)
    out_path = os.path.join(OUT_BASE, cell, f"perseus_{cell}_posterior_bands_v3.csv")
    bands.to_csv(out_path, index=False)
    print(f"    -> wrote {out_path}")

    # Brief summary at year 50 + 95 for NE calibrated
    for var in bands.VARIANT.unique():
        for y in [50, 95]:
            sub = bands[(bands.VARIANT == var) & (bands.PROJ_YEAR == y)]
            if len(sub) > 0:
                med = sub.loc[sub.STAT == "median", "VALUE"].iloc[0]
                q10 = sub.loc[sub.STAT == "q10", "VALUE"].iloc[0]
                q90 = sub.loc[sub.STAT == "q90", "VALUE"].iloc[0]
                n = sub.loc[sub.STAT == "median", "N_DRAWS"].iloc[0]
                print(f"      {var} y{y}: median={med:.0f}, q10={q10:.0f}, q90={q90:.0f}, n_draws={n}")

    return bands


def main():
    print("FVS Posterior v3 Aggregator")
    print("=" * 60)
    all_bands = []
    for cell in cells:
        b = aggregate_cell(cell)
        if b is not None:
            all_bands.append(b)
    if all_bands:
        combined = pd.concat(all_bands, ignore_index=True)
        combined_path = os.path.join(OUT_BASE, "perseus_posterior_bands_v3_all.csv")
        combined.to_csv(combined_path, index=False)
        print(f"\nCombined posterior bands: {combined_path}")
        print(f"  Cells: {combined.SCENARIO.nunique()}, variants: {combined.VARIANT.nunique()}")


if __name__ == "__main__":
    main()
