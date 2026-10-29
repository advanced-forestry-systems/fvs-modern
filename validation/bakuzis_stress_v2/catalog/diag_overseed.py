#!/usr/bin/env python3
"""Diagnose the eastern over-seeding (engine year-0 BA above BALIVE) on seed-gate failures.
For each failing stand: BALIVE, CONDPROP_UNADJ, tree-sum BA for the condition (per plot acre and per
condition acre), number of conditions on the plot, share of saplings (DIA < 5), microplot vs subplot
expansion, and the ratio the engine reported. Public FIA tables only."""
import sys, os, pandas as pd, numpy as np
seed = pd.read_csv(sys.argv[1], dtype={"plt_cn": str})
assign = pd.read_csv(sys.argv[2], dtype={"PLT_CN": str})
fia = os.path.expanduser("~/fia_data")
FIPS2ST = {1: "AL", 5: "AR", 12: "FL", 13: "GA", 21: "KY", 22: "LA", 28: "MS", 37: "NC", 40: "OK", 45: "SC", 47: "TN", 48: "TX", 51: "VA",
           17: "IL", 18: "IN", 19: "IA", 29: "MO", 39: "OH", 20: "KS", 31: "NE", 26: "MI", 27: "MN", 55: "WI", 23: "ME", 33: "NH", 50: "VT", 25: "MA", 9: "CT", 36: "NY", 42: "PA", 44: "RI", 24: "MD", 10: "DE", 54: "WV", 34: "NJ"}
f = seed[(seed.status == "SEED_FAIL") & (seed.config == "default") & (seed.arm == "native")].copy()
f["ratio"] = f.ba0_ft2ac / f.balive
out = []
for v in sys.argv[3:]:
    fv = f[f.variant == v].head(40)
    # recover state via the assignment table: plot CN was recovered from STATECD/COUNTYCD/PLOT/INVYR1 -> use seed plt_cn in PLOT tables
    states = sorted(assign[assign.variant == v].STATECD.unique())
    cond = pd.concat([pd.read_csv(f"{fia}/{FIPS2ST[s]}_COND.csv", usecols=["PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "BALIVE", "STDAGE"], low_memory=False) for s in states if s in FIPS2ST], ignore_index=True)
    cond["PLT_CN"] = cond.PLT_CN.astype(str)
    cns = set(fv.plt_cn)
    plot = pd.concat([pd.read_csv(f"{fia}/{FIPS2ST[s]}_PLOT.csv", usecols=["CN", "STATECD", "COUNTYCD", "PLOT", "INVYR"], low_memory=False) for s in states if s in FIPS2ST], ignore_index=True)
    plot["CN"] = plot.CN.astype(str); plot = plot[plot.CN.isin(cns)]
    key = assign[assign.variant == v].merge(plot, left_on=["STATECD", "COUNTYCD", "PLOT", "INVYR1"], right_on=["STATECD", "COUNTYCD", "PLOT", "INVYR"])
    condid_of = dict(zip(key.CN, key.CONDID))
    tree = pd.concat([pd.read_csv(f"{fia}/{FIPS2ST[s]}_TREE.csv", usecols=["PLT_CN", "CONDID", "STATUSCD", "DIA", "TPA_UNADJ", "SUBP"], low_memory=False).assign(PLT_CN=lambda d: d.PLT_CN.astype(str)).query("PLT_CN in @cns") for s in states if s in FIPS2ST], ignore_index=True)
    for r in fv.itertuples():
        cid = int(condid_of.get(r.plt_cn, 1))
        c = cond[cond.PLT_CN == r.plt_cn]
        t = tree[(tree.PLT_CN == r.plt_cn) & (tree.STATUSCD == 1)]
        tc = t[t.CONDID == cid]
        ba_tree = float((tc.DIA ** 2 * 0.005454154 * tc.TPA_UNADJ).sum())
        ba_sap = float((tc[tc.DIA < 5].DIA ** 2 * 0.005454154 * tc[tc.DIA < 5].TPA_UNADJ).sum())
        out.append(dict(variant=v, plt_cn=r.plt_cn, condid=cid, balive_cond=float(c[c.CONDID == cid].BALIVE.iloc[0]) if (c.CONDID == cid).any() else np.nan, balive=r.balive, ba0=r.ba0_ft2ac, ratio=round(r.ratio, 2), condprop=r.condprop, n_cond=len(c), n_forest_cond=int((c.COND_STATUS_CD == 1).sum()),
                        ba_tree_plotacre=round(ba_tree, 1), ba_tree_condacre=round(ba_tree / max(r.condprop, 1e-6), 1), sap_share=round(ba_sap / max(ba_tree, 1e-9), 2), n_trees=len(tc), max_tpa=float(tc.TPA_UNADJ.max()) if len(tc) else np.nan))
o = pd.DataFrame(out); pd.set_option("display.width", 250); print(o.to_string()); o.to_csv("diag_overseed.csv", index=False)
print(o.groupby("variant")[["ratio", "condprop", "n_cond", "sap_share"]].median())
