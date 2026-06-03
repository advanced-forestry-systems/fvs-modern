#!/usr/bin/env python3
"""treeinit_fix_tpa.py -- repair the FVS TreeInit expansion factor.

ROOT CAUSE (diagnosed 2026-06-03): the DataMart <ST>_FVS_TREEINIT_PLOT.csv files
carry a TREE_COUNT (per-acre expansion) that is ~1/6 of FIA's true TPA_UNADJ for
overstory trees. Maine per-plot basal area came out 34.8 ft2/ac vs the raw FIA
TREE table's 103.5 (median ratio 6.5x), so the whole CONUS campaign biomass was
~6x too light and unpublishable. The join (STAND_CN) and the engine are fine.

FIX: each treeinit row carries TREE_CN (the FIA tree CN). Join it to the raw
<ST>_TREE.csv (CN, TPA_UNADJ, STATUSCD) and overwrite TREE_COUNT with the
authoritative TPA_UNADJ for matched live trees. Everything else (DIAMETER, HT,
CRRATIO, SPECIES, PLOT_ID, STAND_CN) is preserved, so run_conus_task_fvstreeinit
needs no change -- just point --treeinit-dir at the fixed output.

Unmatched rows (older panels not in the current TREE.csv snapshot) keep their
original TREE_COUNT (conservative; leaves them light rather than dropping them).

Usage:
  python3 treeinit_fix_tpa.py --in-dir FIA_fresh/treeinit \
      --tree-dir FIA --out-dir FIA_fresh/treeinit_fixed [--states ME,GA,...]
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd

BA = 0.005454  # ft2 per inch^2 (basal-area constant), for the QA report only


def fix_state(st, in_dir, tree_dir, out_dir):
    tin = os.path.join(in_dir, f"{st}_FVS_TREEINIT_PLOT.csv")
    rtp = os.path.join(tree_dir, f"{st}_TREE.csv")
    if not os.path.exists(tin):
        return st, "no_treeinit", None
    ti = pd.read_csv(tin, dtype={"TREE_CN": str}, low_memory=False)
    if "TREE_CN" not in ti or "TREE_COUNT" not in ti:
        return st, "bad_schema", None
    ti["TREE_COUNT"] = pd.to_numeric(ti["TREE_COUNT"], errors="coerce")
    ti["DIAMETER"] = pd.to_numeric(ti["DIAMETER"], errors="coerce")
    ba_old = float((ti["TREE_COUNT"] * BA * ti["DIAMETER"] ** 2)
                   .groupby(ti["STAND_CN"]).sum().mean())
    matched = np.nan
    if os.path.exists(rtp):
        rc = pd.read_csv(rtp, nrows=0).columns.tolist()
        use = [c for c in ["CN", "TPA_UNADJ", "STATUSCD"] if c in rc]
        if "CN" in use and "TPA_UNADJ" in use:
            rt = pd.read_csv(rtp, usecols=use, dtype={"CN": str}, low_memory=False)
            if "STATUSCD" in rt:                       # live trees only
                rt = rt[rt["STATUSCD"] == 1]
            rt = rt.rename(columns={"CN": "TREE_CN"})[["TREE_CN", "TPA_UNADJ"]]
            rt = rt.dropna(subset=["TREE_CN"]).drop_duplicates("TREE_CN")
            ti = ti.merge(rt, on="TREE_CN", how="left")
            hit = ti["TPA_UNADJ"].notna()
            matched = float(hit.mean())
            ti.loc[hit, "TREE_COUNT"] = ti.loc[hit, "TPA_UNADJ"]
            ti = ti.drop(columns=["TPA_UNADJ"])
    ba_new = float((ti["TREE_COUNT"] * BA * ti["DIAMETER"] ** 2)
                   .groupby(ti["STAND_CN"]).sum().mean())
    os.makedirs(out_dir, exist_ok=True)
    ti.to_csv(os.path.join(out_dir, f"{st}_FVS_TREEINIT_PLOT.csv"), index=False)
    return st, "ok", (ba_old, ba_new, matched)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--tree-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--states", default="")
    a = ap.parse_args()
    if a.states:
        states = [s.strip().upper() for s in a.states.split(",")]
    else:
        states = sorted(os.path.basename(f).split("_")[0]
                        for f in glob.glob(os.path.join(a.in_dir,
                                                        "*_FVS_TREEINIT_PLOT.csv")))
    print(f"{len(states)} states")
    for st in states:
        s, status, q = fix_state(st, a.in_dir, a.tree_dir, a.out_dir)
        if q:
            print(f"  {s}: BA {q[0]:.1f} -> {q[1]:.1f} ft2/ac "
                  f"(TPA matched {100*q[2]:.0f}%)" if q[2] == q[2]
                  else f"  {s}: BA {q[0]:.1f} -> {q[1]:.1f} (no raw TREE)")
        else:
            print(f"  {s}: {status}")


if __name__ == "__main__":
    main()
