#!/usr/bin/env python3
"""ME TreeMap pilot: extract donor trees (from the TreeMap tree table) and donor
standinit locations, for the Maine donor plots from Stage A.

The earlier attempt drew tree lists from FIA_fresh/treeinit, which does not cover
the TreeMap-2022 donor plot CNs (-> near-zero biomass). The correct source for a
TreeMap projection is TreeMap's own tree table, which ships the imputed tree list
per plot. This filters both to the ME donor set.

Outputs:
  me_donor_trees.csv      PLT_CN, SPCD, DIA, HT, CR, TPA_UNADJ, SUBP, TREE (live)
  me_donor_standinit.csv  the ENTIRE-standinit rows for the donor plots (location)
"""
from __future__ import annotations
import argparse
import pandas as pd

TREE_TABLE = ("/fs/scratch/PUOM0008/crsfaaron/reference_rasters/TREEMAP/TM2022/"
              "TreeMap2022_CONUS_Tree_Table.csv")
STANDINIT = "/fs/scratch/PUOM0008/crsfaaron/FIA/ENTIRE_FVS_STANDINIT_PLOT.csv"
SI_KEEP = ["STAND_CN", "VARIANT", "STATE", "INV_YEAR", "LATITUDE", "LONGITUDE",
           "ELEVFT", "SLOPE", "ASPECT", "AGE"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", required=True)
    ap.add_argument("--out-trees", required=True)
    ap.add_argument("--out-standinit", required=True)
    a = ap.parse_args()

    dn = pd.read_csv(a.donors, dtype={"PLT_CN": str})
    want = set(dn["PLT_CN"].dropna().astype(str))
    print(f"{len(want)} donor PLT_CNs", flush=True)

    # ---- trees from the TreeMap tree table (live only) ----
    cols = ["PLT_CN", "STATUSCD", "TPA_UNADJ", "SPCD", "DIA", "HT", "CR",
            "SUBP", "TREE"]
    chunks = []
    for ch in pd.read_csv(TREE_TABLE, usecols=cols, dtype={"PLT_CN": str},
                          low_memory=False, chunksize=500000):
        sel = ch[(ch["STATUSCD"] == 1) & (ch["PLT_CN"].isin(want))]
        if len(sel):
            chunks.append(sel)
    tr = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=cols)
    tr.to_csv(a.out_trees, index=False)
    print(f"trees: {len(tr)} rows, {tr['PLT_CN'].nunique()} plots with trees "
          f"-> {a.out_trees}", flush=True)

    # ---- standinit (locations) for donor plots ----
    sc = []
    for ch in pd.read_csv(STANDINIT, usecols=lambda c: c in SI_KEEP,
                          low_memory=False, chunksize=200000):
        ch["STAND_CN"] = ch["STAND_CN"].apply(
            lambda x: str(int(float(x))) if pd.notna(x) else "")
        sel = ch[ch["STAND_CN"].isin(want)]
        if len(sel):
            sc.append(sel)
    si = pd.concat(sc, ignore_index=True) if sc else pd.DataFrame(columns=SI_KEEP)
    si.to_csv(a.out_standinit, index=False)
    print(f"standinit: {si['STAND_CN'].nunique()} donor plots -> {a.out_standinit}",
          flush=True)


if __name__ == "__main__":
    main()
