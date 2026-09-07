"""Tree-level sanity gates for the balfix rerun.

GATE 2  cycle-0 invariance: every treelist row at Year == MEASYEAR1 must be bit-identical
        to the bin-fixed_20260804 reference run, for every arm. Cycle 0 involves no
        growth so it cannot be touched by a PTBALT repair.
GATE 3  all three Greg arms show nonzero diameter growth and a negative growth-to-BAL
        relationship (BAL quartile mean DG strictly monotone decreasing).
GATE 4b native-arm bit-identity at tree level across the full treelist.
"""
import glob, os, sys, json, numpy as np, pandas as pd

BAL = "/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804/out_full"
FIX = "/fs/scratch/PUOM0008/crsfaaron/finalrun_20260804/out_full"
OUT = "/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805"
ARMS = ["fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate"]
GREG = ["organon", "conus_spdep", "conus_climate"]
NATIVE = ["fvs_base", "fvs_regional"]
UC = ["PLOT", "model", "Year", "TreeId", "TreeIndex", "SpeciesFIA", "TPA", "DBH",
      "DG", "Ht", "HtG", "PctCr", "PtBAL", "MCuFt", "MEASYEAR1", "MEASYEAR2"]
NUM = ["TPA", "DBH", "DG", "Ht", "HtG", "PctCr", "PtBAL", "MCuFt"]
os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 240)
G = {}


def load(d, tag):
    fs = sorted(glob.glob(d + "/treelist_fin_armA_s*.csv"))
    t = pd.concat([pd.read_csv(f, usecols=UC, low_memory=False) for f in fs],
                  ignore_index=True)
    print("[%s] treelist shards=%d rows=%d plots=%d" % (tag, len(fs), len(t), t.PLOT.nunique()))
    sys.stdout.flush()
    return t


tb = load(BAL, "balfix")
tf = load(FIX, "binfixed")

KEY = ["PLOT", "model", "Year", "TreeId"]

print("\n" + "=" * 78)
print("GATE 2. CYCLE-0 INVARIANCE (Year == MEASYEAR1)")
print("=" * 78)
c0b = tb.loc[tb.Year == tb.MEASYEAR1].sort_values(KEY).reset_index(drop=True)
c0f = tf.loc[tf.Year == tf.MEASYEAR1].sort_values(KEY).reset_index(drop=True)
print("cycle-0 rows  balfix=%d  binfixed=%d" % (len(c0b), len(c0f)))
same_shape = (len(c0b) == len(c0f))
keys_same = same_shape and c0b[KEY].equals(c0f[KEY])
print("row counts equal: %s   keys identical: %s" % (same_shape, keys_same))
res0 = {}
if keys_same:
    for arm in ARMS:
        ma = c0b.model == arm
        mx = 0.0; worst = ""
        for c in NUM:
            dd = np.nanmax(np.abs(c0b.loc[ma, c].values.astype(float)
                                  - c0f.loc[ma, c].values.astype(float)))
            if not np.isfinite(dd):
                dd = 0.0
            if dd > mx:
                mx = dd; worst = c
        res0[arm] = mx
        print("  %-15s cycle-0 rows=%d  maxabs=%.6g %s"
              % (arm, int(ma.sum()), mx, ("(" + worst + ")") if worst else ""))
    G["gate2_cycle0_invariant"] = all(v == 0.0 for v in res0.values())
else:
    G["gate2_cycle0_invariant"] = False
print("GATE 2 PASS: %s" % G["gate2_cycle0_invariant"])

print("\n" + "=" * 78)
print("GATE 4b. NATIVE-ARM TREE-LEVEL BIT-IDENTITY, FULL TREELIST")
print("=" * 78)
ab = tb.sort_values(KEY).reset_index(drop=True)
af = tf.sort_values(KEY).reset_index(drop=True)
resall = {}
if len(ab) == len(af) and ab[KEY].equals(af[KEY]):
    for arm in ARMS:
        ma = ab.model == arm
        mx = 0.0; worst = ""
        for c in NUM:
            dd = np.nanmax(np.abs(ab.loc[ma, c].values.astype(float)
                                  - af.loc[ma, c].values.astype(float)))
            if not np.isfinite(dd):
                dd = 0.0
            if dd > mx:
                mx = dd; worst = c
        resall[arm] = mx
        print("  %-15s rows=%d  maxabs=%.6g %s"
              % (arm, int(ma.sum()), mx, ("(" + worst + ")") if worst else ""))
    G["gate4b_native_identical"] = all(resall[a] == 0.0 for a in NATIVE)
    G["gate4b_greg_moved"] = all(resall[a] > 0.0 for a in GREG)
else:
    print("treelist key sets differ: balfix rows=%d binfixed rows=%d" % (len(ab), len(af)))
    G["gate4b_native_identical"] = False
    G["gate4b_greg_moved"] = False
print("GATE 4b native bit-identical: %s   Greg moved: %s"
      % (G["gate4b_native_identical"], G["gate4b_greg_moved"]))
if not G["gate4b_native_identical"]:
    print("\n*** STOP. NATIVE ARM MOVED AT TREE LEVEL. LEAKAGE. ***")

print("\n" + "=" * 78)
print("GATE 3. GREG ARMS: NONZERO DIAMETER GROWTH AND NEGATIVE GROWTH-TO-BAL")
print("=" * 78)
gr = tb.loc[(tb.Year > tb.MEASYEAR1) & tb.DG.notna()].copy()
print("growth rows = %d" % len(gr))
rows = []
mono = {}
for arm in ARMS:
    g = gr.loc[gr.model == arm]
    if not len(g):
        continue
    dgm = float(g.DG.mean()); nz = float((g.DG > 0).mean())
    q = pd.qcut(g.PtBAL.rank(method="first"), 4, labels=[1, 2, 3, 4])
    qm = g.groupby(q, observed=True).DG.mean()
    qb = g.groupby(q, observed=True).PtBAL.mean()
    corr = float(np.corrcoef(g.PtBAL.values.astype(float),
                             g.DG.values.astype(float))[0, 1])
    sp = float(pd.Series(g.PtBAL.values).corr(pd.Series(g.DG.values), method="spearman"))
    vals = [float(v) for v in qm.values]
    isdec = all(vals[i] > vals[i + 1] for i in range(3))
    mono[arm] = isdec
    print("  %-15s meanDG=%.5f  frac DG>0=%.4f  r=%+.4f  rho=%+.4f  monotone_dec=%s"
          % (arm, dgm, nz, corr, sp, isdec))
    print("      BAL quartile mean DG : " + ", ".join("%.4f" % v for v in vals))
    print("      BAL quartile mean BAL: " + ", ".join("%.1f" % v for v in qb.values))
    rows.append(dict(ARM=arm, meanDG=dgm, frac_pos=nz, pearson=corr, spearman=sp,
                     q1=vals[0], q2=vals[1], q3=vals[2], q4=vals[3], monotone=isdec))
pd.DataFrame(rows).to_csv(OUT + "/greg_dg_bal.csv", index=False)
G["gate3_greg_dg_nonzero"] = all(
    float(gr.loc[gr.model == a].DG.mean()) > 0 for a in GREG)
G["gate3_greg_bal_negative"] = all(mono.get(a, False) for a in GREG)
print("GATE 3 nonzero DG on all Greg arms : %s" % G["gate3_greg_dg_nonzero"])
print("GATE 3 monotone decreasing DG by BAL quartile on all Greg arms : %s"
      % G["gate3_greg_bal_negative"])

with open(OUT + "/sanity_treelist.json", "w") as f:
    json.dump({k: bool(v) for k, v in G.items()}, f, indent=1)
print("\nSANITY (tree side):", {k: bool(v) for k, v in G.items()})
