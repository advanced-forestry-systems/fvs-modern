#!/usr/bin/env python3
"""ME TreeMap pilot, Stage B (v2): run FVS on the Maine donor plots using
TreeMap's own tree lists.

Builds the FVS treeinit from the TreeMap-2022 tree table (me_donor_trees.csv,
the imputed donor tree list -- the correct, self-consistent source for a TreeMap
projection) and the standinit location from me_donor_standinit.csv, then projects
100 yr in the requested configs. Output me_donor_trajectories[...].csv:
    PLT_CN, STATE, VARIANT, CONFIG, PROJ_YEAR, AGB_TONS_AC

Gompit arm: set FVS_LIB_DIR=fvs_gompit/lib + FVS_GOMPIT=1 and --configs default,
then relabel CONFIG=gompit downstream.
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


def treeinit_from_treemap(trees: pd.DataFrame, sid: str) -> pd.DataFrame:
    """TreeMap tree table rows -> fvs_treeinit schema."""
    recs = []
    for i, t in enumerate(trees.itertuples(index=False)):
        dia = getattr(t, "DIA", np.nan)
        if pd.isna(dia) or float(dia) < 1.0:
            continue
        ht = getattr(t, "HT", 0)
        cr = getattr(t, "CR", 0)
        recs.append({
            "stand_id": sid,
            "plot_id": int(float(getattr(t, "SUBP", 1) or 1)),
            "tree_id": i + 1,
            "tree_count": float(getattr(t, "TPA_UNADJ", 1.0) or 1.0),
            "species": int(float(getattr(t, "SPCD", 0) or 0)),
            "diameter": round(float(dia), 1),
            "ht": round(float(ht), 0) if pd.notna(ht) and float(ht) > 0 else 0,
            "crratio": int(float(cr)) if pd.notna(cr) and float(cr) > 0 else 0,
        })
    return pd.DataFrame(recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--standinit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--configs", default="default,calibrated")
    a = ap.parse_args()
    cfgs = [c.strip() for c in a.configs.split(",")]
    arm = [(None if c == "default" else c, c) for c in cfgs]

    tr = pd.read_csv(a.trees, dtype={"PLT_CN": str})
    by_plot = {k: v for k, v in tr.groupby("PLT_CN")}
    si = pd.read_csv(a.standinit, low_memory=False)
    si["STAND_CN"] = si["STAND_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else "")
    print(f"{len(by_plot)} plots with TreeMap trees, {len(si)} standinit rows",
          flush=True)

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    rows = []
    n = 0
    for _, stand in si.iterrows():
        cn = stand["STAND_CN"]
        trees = by_plot.get(cn)
        if trees is None or trees.empty:
            continue
        variant = str(stand.get("VARIANT", "ne")).lower()
        try:
            state = FIPS[int(float(stand.get("STATE")))]
        except (KeyError, ValueError, TypeError):
            state = "ME"
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
            tdf = treeinit_from_treemap(trees, sid)
        except Exception:
            continue
        if tdf.empty:
            continue
        n += 1
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
        if n % 500 == 0:
            print(f"  {n} plots projected", flush=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"projected {n} plots, wrote {len(rows)} rows -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
