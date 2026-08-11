#!/usr/bin/env python3
"""treeinit_fix_v2.py -- treeinit repair v2: TPA_UNADJ expansion + complete heights.

v1 (treeinit_fix_tpa.py) fixed the 6.5x under-expanded TREE_COUNT. This v2 also
completes tree HEIGHTS, because FVS reports the PROJ_YEAR 0 (2025) inventory
treelist before it imputes missing heights, and NSBE biomass needs height -- so
the 2025 point was understated ~50-90% (the height-fill artifact we dodged by
anchoring at 2030). Filling heights up front makes 2025 whole and lets us anchor
the dashboard at the true inventory year.

Per state, joining each treeinit row to the raw FIA TREE table by TREE_CN:
  TREE_COUNT <- TPA_UNADJ            (expansion fix, v1)
  HT         <- TREE.HT where the treeinit HT is missing (FIA modeled height,
                ~85% of live trees)
Residual missing heights (~15%) are imputed from a per-species log-log H-D fit
(HT = exp(a + b*log(DIA))) on the trees that have height, with a pooled fallback;
predictions clamped to [5, 320] ft.

Output: <out>/<ST>_FVS_TREEINIT_PLOT.csv  (engine-ready, run unchanged).
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd

BA = 0.005454
HT_MIN, HT_MAX = 5.0, 320.0


def hd_impute(df):
    """fill HT<=0 from per-species log-log H-D fit; pooled fallback."""
    have = df["HT"] > 0
    need = ~have & (df["DIAMETER"] > 0)
    if need.sum() == 0:
        return df
    ld = np.log(df["DIAMETER"].clip(lower=1.0))
    lh = np.log(df["HT"].clip(lower=1.0))
    # pooled fit
    h = have & np.isfinite(ld) & np.isfinite(lh)
    if h.sum() >= 20:
        pb, pa = np.polyfit(ld[h], lh[h], 1)
    else:
        pb, pa = 0.5, 2.5
    pred = pd.Series(np.exp(pa + pb * ld), index=df.index)
    # per-species refinement where enough data
    for sp, g in df[have].groupby("SPECIES"):
        if len(g) >= 30:
            x, y = np.log(g["DIAMETER"].clip(lower=1.0)), np.log(g["HT"])
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 30:
                b, a = np.polyfit(x[ok], y[ok], 1)
                m = need & (df["SPECIES"] == sp)
                pred.loc[m] = np.exp(a + b * ld[m])
    df.loc[need, "HT"] = pred.loc[need].clip(HT_MIN, HT_MAX).round(0)
    return df


def fix_state(st, in_dir, tree_dir, out_dir):
    tin = os.path.join(in_dir, f"{st}_FVS_TREEINIT_PLOT.csv")
    rtp = os.path.join(tree_dir, f"{st}_TREE.csv")
    if not os.path.exists(tin):
        return st, "no_treeinit", None
    ti = pd.read_csv(tin, dtype={"TREE_CN": str}, low_memory=False)
    for c in ("TREE_COUNT", "DIAMETER", "HT"):
        if c in ti:
            ti[c] = pd.to_numeric(ti[c], errors="coerce")
    ht0 = float((ti["HT"] > 0).mean())
    matched = np.nan
    if os.path.exists(rtp):
        rc = pd.read_csv(rtp, nrows=0).columns.tolist()
        use = [c for c in ["CN", "TPA_UNADJ", "HT", "STATUSCD"] if c in rc]
        if "CN" in use and "TPA_UNADJ" in use:
            rt = pd.read_csv(rtp, usecols=use, dtype={"CN": str}, low_memory=False)
            if "STATUSCD" in rt:
                rt = rt[rt["STATUSCD"] == 1]
            rt = rt.rename(columns={"CN": "TREE_CN", "HT": "HT_FIA"})
            keep = ["TREE_CN", "TPA_UNADJ"] + (["HT_FIA"] if "HT_FIA" in rt else [])
            rt = rt[keep].dropna(subset=["TREE_CN"]).drop_duplicates("TREE_CN")
            ti = ti.merge(rt, on="TREE_CN", how="left")
            hit = ti["TPA_UNADJ"].notna()
            matched = float(hit.mean())
            ti.loc[hit, "TREE_COUNT"] = ti.loc[hit, "TPA_UNADJ"]
            if "HT_FIA" in ti:                      # fill missing heights from FIA
                fillh = (ti["HT"].isna() | (ti["HT"] <= 0)) & (ti["HT_FIA"] > 0)
                ti.loc[fillh, "HT"] = ti.loc[fillh, "HT_FIA"]
            ti = ti.drop(columns=[c for c in ("TPA_UNADJ", "HT_FIA") if c in ti])
    ti["HT"] = ti["HT"].fillna(0)
    ti = hd_impute(ti)                              # impute the residual
    ht1 = float((ti["HT"] > 0).mean())
    os.makedirs(out_dir, exist_ok=True)
    ti.to_csv(os.path.join(out_dir, f"{st}_FVS_TREEINIT_PLOT.csv"), index=False)
    return st, "ok", (ht0, ht1, matched)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--tree-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--states", default="")
    a = ap.parse_args()
    states = ([s.strip().upper() for s in a.states.split(",")] if a.states
              else sorted(os.path.basename(f).split("_")[0]
                          for f in glob.glob(os.path.join(a.in_dir,
                                                          "*_FVS_TREEINIT_PLOT.csv"))))
    print(f"{len(states)} states")
    for st in states:
        s, status, q = fix_state(st, a.in_dir, a.tree_dir, a.out_dir)
        if q:
            print(f"  {s}: HT coverage {100*q[0]:.0f}% -> {100*q[1]:.0f}% "
                  f"(TPA matched {100*q[2]:.0f}%)" if q[2] == q[2]
                  else f"  {s}: HT {100*q[0]:.0f}%->{100*q[1]:.0f}% (no raw TREE)")
        else:
            print(f"  {s}: {status}")


if __name__ == "__main__":
    main()
