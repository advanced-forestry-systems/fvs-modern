#!/usr/bin/env python3
"""
CONUS-wide 100-year no-harvest FVS projection (default vs calibrated), all variants.

Sources stands from the pre-built CONUS FVS stand-init
(ENTIRE_FVS_STANDINIT_COND.csv, which carries the FVS VARIANT per stand) and
trees from the per-state FIA TREE CSVs. Reuses the perseus engine internals
(run_fvs_projection, build_fvs_standinit, build_fvs_treeinit, compute_plot_agb,
NSBECalculator) so the calibrated runs pick up the calibration_multipliers block
via FvsConfigLoader.generate_keywords().

Usage (per-variant SLURM array task):
  python conus_100yr_projection.py --variant CS --sample 300 --output-dir out/
  python conus_100yr_projection.py --variant SN --batch-id 3 --batch-size 5000 --output-dir out/

Designed for SLURM array jobs: one task per variant (stress-test sample) or per
(variant, batch) for the full CONUS run.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
import perseus_100yr_projection as P  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("conus")

FIA_DIR = os.environ.get("FVS_FIA_DATA_DIR", "/fs/scratch/PUOM0008/crsfaaron/FIA")
# PLOT-level stand-init: STAND_CN equals the FIA tree PLT_CN (the COND-level file
# is keyed by condition CN, which does not join to the per-state tree files).
STANDINIT = os.path.join(FIA_DIR, "ENTIRE_FVS_STANDINIT_PLOT.csv")

FIPS = {1: "AL", 2: "AK", 4: "AZ", 5: "AR", 6: "CA", 8: "CO", 9: "CT", 10: "DE",
        12: "FL", 13: "GA", 16: "ID", 17: "IL", 18: "IN", 19: "IA", 20: "KS",
        21: "KY", 22: "LA", 23: "ME", 24: "MD", 25: "MA", 26: "MI", 27: "MN",
        28: "MS", 29: "MO", 30: "MT", 31: "NE", 32: "NV", 33: "NH", 34: "NJ",
        35: "NM", 36: "NY", 37: "NC", 38: "ND", 39: "OH", 40: "OK", 41: "OR",
        42: "PA", 44: "RI", 45: "SC", 46: "SD", 47: "TN", 48: "TX", 49: "UT",
        50: "VT", 51: "VA", 53: "WA", 54: "WV", 55: "WI", 56: "WY"}

KEEP = ["STAND_CN", "STAND_ID", "VARIANT", "STATE", "INV_YEAR",
        "LATITUDE", "LONGITUDE", "ELEVFT", "SLOPE", "ASPECT", "AGE"]


def load_variant_stands(variant, sample=None, batch_id=0, batch_size=None):
    """Read ENTIRE standinit (chunked), filter to one FVS variant, sample/batch."""
    want = variant.upper()
    chunks = []
    for ch in pd.read_csv(STANDINIT, usecols=lambda c: c in KEEP,
                          low_memory=False, chunksize=200000):
        sel = ch[ch["VARIANT"].astype(str).str.upper() == want]
        if len(sel):
            chunks.append(sel)
    si = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=KEEP)
    if sample and len(si) > sample:
        si = si.sample(n=sample, random_state=42).reset_index(drop=True)
    elif batch_size:
        si = si.iloc[batch_id * batch_size:(batch_id + 1) * batch_size].reset_index(drop=True)
    return si


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", required=True, help="FVS variant code (e.g. CS, SN, NE)")
    ap.add_argument("--sample", type=int, default=None,
                    help="random stand sample size (stress test)")
    ap.add_argument("--batch-id", type=int, default=0)
    ap.add_argument("--batch-size", type=int, default=None,
                    help="full-run batch size (mutually exclusive with --sample)")
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.output_dir, exist_ok=True)

    si = load_variant_stands(a.variant, a.sample, a.batch_id, a.batch_size)
    log.info(f"variant {a.variant}: {len(si)} stands selected")
    if si.empty:
        log.warning("no stands; exiting")
        return

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    rows = []
    n_done = 0
    for state_fips, grp in si.groupby("STATE"):
        try:
            state = FIPS[int(float(state_fips))]
        except (KeyError, ValueError, TypeError):
            log.warning(f"unknown STATE={state_fips}; skipping {len(grp)} stands")
            continue
        plt_cns = [str(int(float(x))) for x in grp["STAND_CN"]]
        try:
            trees = P.load_fia_trees_for_plots(plt_cns, FIA_DIR, state)
        except Exception as e:
            log.warning(f"state {state}: tree load failed: {e}")
            continue
        if trees.empty or "PLT_CN" not in trees.columns:
            log.info(f"state {state}: no matching trees for {len(grp)} stands")
            continue
        tcn = trees["PLT_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
        for _, stand in grp.iterrows():
            cn = str(int(float(stand["STAND_CN"])))
            st = trees.loc[tcn == cn]
            if st.empty:
                continue
            sid = f"S{cn}"
            inv_year = int(float(stand.get("INV_YEAR") or 2010))
            plot_data = {
                "INVYR": inv_year,
                "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
                "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
                "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50,
            }
            try:
                sdf = P.build_fvs_standinit(plot_data, sid, a.variant.lower())
                tdf = P.build_fvs_treeinit(st, sid)
            except Exception as e:
                log.warning(f"{cn}: build failed: {e}")
                continue
            for cfg, label in ((None, "default"), ("calibrated", "calibrated")):
                try:
                    fr = P.run_fvs_projection(sdf, tdf, sid, a.variant.lower(),
                                              config_version=cfg,
                                              num_cycles=20, cycle_length=5)
                    for cy, tl in sorted(fr["treelists"].items()):
                        py = cy - inv_year
                        if py < 0:
                            continue
                        agb = P.compute_plot_agb(tl, nsbe)
                        rows.append({"STAND_CN": cn, "STATE": state, "YEAR": cy,
                                     "PROJ_YEAR": py, "VARIANT": a.variant.upper(),
                                     "CONFIG": label, "AGB_TONS_AC": round(float(agb), 4)})
                except Exception as e:
                    log.error(f"{cn}/{label}: FVS failed: {e}")
            n_done += 1
            if n_done % 50 == 0:
                log.info(f"  {n_done} stands projected")

    tag = f"s{a.sample}" if a.sample else f"b{a.batch_id}"
    out = os.path.join(a.output_dir, f"conus_{a.variant.lower()}_{tag}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    log.info(f"variant {a.variant}: {n_done} stands, {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
