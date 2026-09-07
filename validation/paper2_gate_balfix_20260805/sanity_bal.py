"""GATE 3, done properly.

The treelist PtBAL column is written as identically 0.0 by the driver, so it cannot be
used. BAL is recomputed here from the cycle-0 tree list: for each (PLOT, model) stand,
BAL for a tree is the sum of basal area per acre of all trees strictly larger in DBH.
Diameter growth is the treelist DG on the projection row, annualized by PERIOD_YR.

Reports, per arm and per binary (balfix vs bin-fixed reference):
  mean annual DG, fraction of trees with DG > 0, Pearson and Spearman correlation of
  annual DG against cycle-0 BAL, and mean annual DG by BAL quartile.
"""
import glob, os, sys, json, numpy as np, pandas as pd

BAL = "/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804/out_full"
FIX = "/fs/scratch/PUOM0008/crsfaaron/finalrun_20260804/out_full"
OUT = "/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805"
ARMS = ["fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate"]
GREG = ["organon", "conus_spdep", "conus_climate"]
UC = ["PLOT", "model", "Year", "TreeId", "TPA", "DBH", "DG", "PERIOD_YR", "MEASYEAR1"]
pd.set_option("display.width", 240)


def prep(d, tag):
    fs = sorted(glob.glob(d + "/treelist_fin_armA_s*.csv"))
    t = pd.concat([pd.read_csv(f, usecols=UC, low_memory=False) for f in fs],
                  ignore_index=True)
    print("[%s] rows=%d" % (tag, len(t))); sys.stdout.flush()
    c0 = t.loc[t.Year == t.MEASYEAR1].copy()
    gr = t.loc[t.Year > t.MEASYEAR1].copy()
    # BAL at cycle 0: sum of BA/ac of all strictly larger trees in the same stand
    c0["BAt"] = 0.005454154 * c0.DBH.astype(float) ** 2 * c0.TPA.astype(float)
    c0 = c0.sort_values(["PLOT", "model", "DBH"], ascending=[True, True, False])
    g = c0.groupby(["PLOT", "model"], sort=False)
    c0["BAL"] = g.BAt.cumsum() - c0.BAt
    c0["BA_STAND"] = g.BAt.transform("sum")
    c0["RELBAL"] = c0.BAL / c0.BA_STAND.replace(0, np.nan)
    k = ["PLOT", "model", "TreeId"]
    m = gr.merge(c0[k + ["BAL", "RELBAL", "BA_STAND", "DBH"]].rename(
        columns={"DBH": "DBH0"}), on=k, how="inner")
    m["DG_ANN"] = m.DG.astype(float) / m.PERIOD_YR.astype(float)
    print("[%s] matched growth rows=%d  mean BAL=%.2f  mean DG_ANN=%.5f"
          % (tag, len(m), m.BAL.mean(), m.DG_ANN.mean()))
    return m


rows = []
res = {}
for tag, d in [("balfix", BAL), ("binfixed", FIX)]:
    m = prep(d, tag)
    print("\n" + "=" * 78)
    print("BINARY: %s" % tag)
    print("=" * 78)
    for a in ARMS:
        g = m.loc[m.model == a]
        q = pd.qcut(g.BAL.rank(method="first"), 4, labels=[1, 2, 3, 4])
        qm = g.groupby(q, observed=True).DG_ANN.mean()
        qb = g.groupby(q, observed=True).BAL.mean()
        v = [float(x) for x in qm.values]
        pr = float(np.corrcoef(g.BAL.values.astype(float),
                               g.DG_ANN.values.astype(float))[0, 1])
        sp = float(pd.Series(g.BAL.values).corr(pd.Series(g.DG_ANN.values),
                                                method="spearman"))
        mono = all(v[i] > v[i + 1] for i in range(3))
        print("  %-15s n=%d meanDGann=%.5f fracDG>0=%.4f r=%+.4f rho=%+.4f monotone_dec=%s"
              % (a, len(g), g.DG_ANN.mean(), (g.DG > 0).mean(), pr, sp, mono))
        print("      BAL quartile mean annual DG : " + ", ".join("%.4f" % x for x in v))
        print("      BAL quartile mean BAL       : " + ", ".join("%.1f" % x for x in qb.values))
        rows.append(dict(BIN=tag, ARM=a, n=len(g), meanDG_ann=float(g.DG_ANN.mean()),
                         frac_pos=float((g.DG > 0).mean()), pearson=pr, spearman=sp,
                         q1=v[0], q2=v[1], q3=v[2], q4=v[3], monotone=mono))
        if tag == "balfix":
            res[a] = dict(mono=mono, pr=pr, meanDG=float(g.DG_ANN.mean()))

r = pd.DataFrame(rows)
r.to_csv(OUT + "/gate3_bal_response.csv", index=False)
G = {}
G["gate3_greg_dg_nonzero"] = all(res[a]["meanDG"] > 0 for a in GREG)
G["gate3_greg_bal_negative_corr"] = all(res[a]["pr"] < 0 for a in GREG)
G["gate3_greg_bal_monotone"] = all(res[a]["mono"] for a in GREG)
G["gate3_PASS"] = G["gate3_greg_dg_nonzero"] and G["gate3_greg_bal_negative_corr"]
print("\nGATE 3 nonzero DG on all Greg arms        : %s" % G["gate3_greg_dg_nonzero"])
print("GATE 3 negative DG-to-BAL correlation      : %s" % G["gate3_greg_bal_negative_corr"])
print("GATE 3 monotone decreasing by BAL quartile : %s" % G["gate3_greg_bal_monotone"])
print("GATE 3 PASS                                : %s" % G["gate3_PASS"])
with open(OUT + "/sanity_bal.json", "w") as f:
    json.dump({k: bool(v) for k, v in G.items()}, f, indent=1)
print("\nFULL TABLE")
print(r.round(5).to_string(index=False))
