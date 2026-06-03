#!/usr/bin/env python3
"""ME TreeMap pilot, Stage B: run FVS on the Maine donor FIA plots.

Reads me_treemap_donors.csv (PLT_CN list from Stage A), finds each plot in the
CONUS standinit (STATE, VARIANT, location), builds standinit + treeinit, and
projects 100 yr in default and calibrated FVS (add gompit via a second pass with
the gompit lib + FVS_GOMPIT=1). Output me_donor_trajectories.csv:
    PLT_CN, STATE, VARIANT, CONFIG, PROJ_YEAR, AGB_TONS_AC

Stage C then paints ME = sum over donors of AGB_density(year) x area_ha (TreeMap
expansion) and compares to FIADB expansion.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P
import run_conus_task_fvstreeinit as RC

FIPS = RC.FIPS
STANDINIT = "/fs/scratch/PUOM0008/crsfaaron/FIA/ENTIRE_FVS_STANDINIT_PLOT.csv"
TI_DIR = "/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit"
KEEP = ["STAND_CN", "VARIANT", "STATE", "INV_YEAR", "LATITUDE", "LONGITUDE",
        "ELEVFT", "SLOPE", "ASPECT", "AGE"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--configs", default="default,calibrated")
    a = ap.parse_args()
    cfgs = [c.strip() for c in a.configs.split(",")]
    arm = [(None if c == "default" else c, c) for c in cfgs]

    dn = pd.read_csv(a.donors, dtype={"PLT_CN": str})
    want = set(dn["PLT_CN"].dropna().astype(str))
    print(f"{len(want)} donor PLT_CNs", flush=True)

    # pull donor stands from the CONUS standinit
    chunks = []
    for ch in pd.read_csv(STANDINIT, usecols=lambda c: c in KEEP,
                          low_memory=False, chunksize=200000):
        ch["STAND_CN"] = ch["STAND_CN"].apply(
            lambda x: str(int(float(x))) if pd.notna(x) else "")
        sel = ch[ch["STAND_CN"].isin(want)]
        if len(sel):
            chunks.append(sel)
    si = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=KEEP)
    print(f"matched {si['STAND_CN'].nunique()} of {len(want)} donors in standinit",
          flush=True)

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    rows = []
    cache = {}
    for state_fips, grp in si.groupby("STATE"):
        try:
            state = FIPS[int(float(state_fips))]
        except (KeyError, ValueError, TypeError):
            continue
        tfile = os.path.join(TI_DIR, f"{state}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(tfile):
            continue
        if state not in cache:
            tt = pd.read_csv(tfile, low_memory=False)
            tt["STAND_CN"] = tt["STAND_CN"].apply(
                lambda x: str(int(float(x))) if pd.notna(x) else "")
            cache[state] = {k: v for k, v in tt.groupby("STAND_CN")}
        by_cn = cache[state]
        for _, stand in grp.iterrows():
            cn = stand["STAND_CN"]
            fvs_rows = by_cn.get(cn)
            if fvs_rows is None or fvs_rows.empty:
                continue
            variant = str(stand.get("VARIANT", "ne")).lower()
            sid = f"S{cn}"
            iy = int(float(stand.get("INV_YEAR") or 2010))
            pdat = {"INVYR": iy, "LAT": stand.get("LATITUDE"),
                    "LON": stand.get("LONGITUDE"),
                    "ELEV": stand.get("ELEVFT") or 500,
                    "SLOPE": stand.get("SLOPE") or 10,
                    "ASPECT": stand.get("ASPECT") or 180,
                    "STDAGE": stand.get("AGE") or 50}
            try:
                sdf = P.build_fvs_standinit(pdat, sid, variant)
                tdf = RC.treeinit_for_stand(fvs_rows, sid)
            except Exception:
                continue
            if tdf.empty:
                continue
            for cfg, label in arm:
                try:
                    fr = P.run_fvs_projection(sdf, tdf, sid, variant,
                                              config_version=cfg,
                                              num_cycles=20, cycle_length=5)
                    for cy, tl in sorted(fr["treelists"].items()):
                        py = cy - iy
                        if py < 0:
                            continue
                        agb = P.compute_plot_agb(tl, nsbe)
                        rows.append({"PLT_CN": cn, "STATE": state,
                                     "VARIANT": variant.upper(), "CONFIG": label,
                                     "PROJ_YEAR": py,
                                     "AGB_TONS_AC": round(float(agb), 4)})
                except Exception:
                    pass
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"wrote {len(rows)} rows -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
