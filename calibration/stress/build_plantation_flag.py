#!/usr/bin/env python3
"""build_plantation_flag.py -- per-plot plantation flag from FIADB COND.

A plot is a plantation if the majority (by CONDPROP_UNADJ) of its forested area
is in conditions with STDORGCD == 1 (clear evidence of artificial regeneration).
This drives the plantation-aware management scenarios (intensive harvest confined
to plantations; see fvs_managed_v2.py).

Output: plt_plantation.csv  (PLT_CN, plantation 0/1, STATE).
National plantation fraction ~10% (GA 28%, MS 27%, OR 20%, ME 2.4%).
"""
import glob
import os

import pandas as pd

COND = os.environ.get("FIA_COND_DIR", "/fs/scratch/PUOM0008/crsfaaron/FIA")
NEED = ["PLT_CN", "COND_STATUS_CD", "STDORGCD", "CONDPROP_UNADJ"]


def main():
    rows, skipped = [], []
    for f in sorted(glob.glob(os.path.join(COND, "*_COND.csv"))):
        st = os.path.basename(f).split("_")[0]
        cols = pd.read_csv(f, nrows=0).columns.tolist()
        if not all(c in cols for c in NEED):
            skipped.append(st)
            continue
        d = pd.read_csv(f, usecols=NEED, dtype={"PLT_CN": str}, low_memory=False)
        d = d[d.COND_STATUS_CD == 1]                      # forested conditions
        d = d.assign(pl_prop=(d.STDORGCD == 1) * d.CONDPROP_UNADJ)
        g = d.groupby("PLT_CN").agg(pl=("pl_prop", "sum"),
                                    tot=("CONDPROP_UNADJ", "sum"))
        g["plantation"] = (g.pl > 0.5 * g.tot).astype(int)
        pf = g.reset_index()[["PLT_CN", "plantation"]]
        pf["STATE"] = st
        rows.append(pf)
    allp = pd.concat(rows, ignore_index=True)
    allp["PLT_CN"] = allp.PLT_CN.str.replace(r"\.0$", "", regex=True)
    allp.to_csv("plt_plantation.csv", index=False)
    print(f"plots {len(allp)} plantation {int(allp.plantation.sum())} "
          f"{100 * allp.plantation.mean():.1f}%  skipped: {skipped}")


if __name__ == "__main__":
    main()
