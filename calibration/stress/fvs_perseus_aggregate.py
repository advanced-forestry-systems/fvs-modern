#!/usr/bin/env python3
"""Aggregate the CONUS FVS campaign output into PERSEUS ycx state-series CSVs.

FVS (default / calibrated / gompit) is a drop-in growth engine for the existing
PERSEUS yield-curve pipeline (perseus_wire/scripts/yield_curve_engine). That
pipeline's merge step (ycx_merge_perseus.py) reads, per state, a
`ycx_<ST>_state_series.csv` with rows

    metric, mgmt, year, value, value_lo, value_hi

in native per-ha density, then expands to state totals via n_plots x A0 (A0
anchored to fia.json carbon). This script produces those CSVs from the campaign
`conus_<variant>_b<batch>.csv` files (columns STAND_CN, STATE, YEAR, PROJ_YEAR,
VARIANT, CONFIG, AGB_TONS_AC).

Mapping:
  * calendar year = START(=2025) + PROJ_YEAR  (align every stand's t0 to the
    PERSEUS baseline year, so all states share a 2025-2125 axis).
  * density: AGB_TONS_AC (short tons/ac) -> Mg/ha  (x 2.241702).
  * metric agb_dry      = AGB density (Mg/ha).
  * metric agc_live_total = AGB density x 0.47 (Mg C/ha).
  * value = across-plot mean; value_lo/hi = 10th/90th percentile (per-plot
    spread within the state-year).
  * mgmt = "reserve (no harvest)" (the campaign is the passive-succession run;
    the managed-harvest bucket comes later via the conus_hcs coupling).

One engine per CONFIG: default->fvs_default, calibrated->fvs_calibrated,
gompit->fvs_gompit. A counts.csv carries n_plots per state (for the merge).

Usage:
  python fvs_perseus_aggregate.py --in-dir out_fvs --out-dir perseus_series \
      --config calibrated --start 2025
"""
from __future__ import annotations
import argparse
import glob
import os

import numpy as np
import pandas as pd

TONS_AC_TO_MG_HA = 2.241702      # short tons/ac -> Mg/ha
C_FRACTION = 0.47                # dry biomass -> carbon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True, help="dir of conus_<v>_b<b>.csv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="calibrated",
                    choices=["default", "calibrated", "gompit"])
    ap.add_argument("--start", type=int, default=2025)
    ap.add_argument("--engine", default=None,
                    help="engine model name; default fvs_<config>")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    engine = a.engine or f"fvs_{a.config}"

    files = sorted(glob.glob(os.path.join(a.in_dir, "conus_*.csv")))
    if not files:
        raise SystemExit(f"no conus_*.csv in {a.in_dir}")
    frames = []
    for f in files:
        try:
            d = pd.read_csv(f, usecols=["STAND_CN", "STATE", "PROJ_YEAR",
                                        "CONFIG", "AGB_TONS_AC"])
        except Exception:
            continue
        frames.append(d[d["CONFIG"] == a.config])
    df = pd.concat(frames, ignore_index=True)
    if df.empty:
        raise SystemExit(f"no rows for config {a.config}")

    df["year"] = a.start + df["PROJ_YEAR"].astype(int)
    df["dens_agb"] = df["AGB_TONS_AC"] * TONS_AC_TO_MG_HA      # Mg/ha
    df["dens_agc"] = df["dens_agb"] * C_FRACTION               # Mg C/ha

    counts = []
    for st, g in df.groupby("STATE"):
        rows = []
        nplots = g["STAND_CN"].nunique()
        counts.append({"state": st, "n_plots": nplots, "engine": engine})
        for metric, col in (("agb_dry", "dens_agb"), ("agc_live_total", "dens_agc")):
            agg = (g.groupby("year")[col]
                     .agg(value="mean",
                          value_lo=lambda x: np.percentile(x, 10),
                          value_hi=lambda x: np.percentile(x, 90))
                     .reset_index())
            for _, r in agg.iterrows():
                rows.append({"metric": metric, "mgmt": "reserve (no harvest)",
                             "year": int(r["year"]),
                             "value": round(float(r["value"]), 4),
                             "value_lo": round(float(r["value_lo"]), 4),
                             "value_hi": round(float(r["value_hi"]), 4)})
        out = pd.DataFrame(rows).sort_values(["metric", "year"])
        out.to_csv(os.path.join(a.out_dir, f"ycx_{st}_state_series.csv"),
                   index=False)
    pd.DataFrame(counts).to_csv(
        os.path.join(a.out_dir, f"fvs_{a.config}_state_counts.csv"), index=False)
    print(f"engine {engine}: {df['STATE'].nunique()} states, "
          f"{df['STAND_CN'].nunique()} plots, "
          f"years {int(df['year'].min())}-{int(df['year'].max())} -> {a.out_dir}")


if __name__ == "__main__":
    main()
