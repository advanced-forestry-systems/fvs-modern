"""Definitive balfix gate on BOTH cohorts.

Arithmetic copied verbatim from finalrun_20260804/gate_final.py (itself verbatim from
paper2_volgate_20260728/gate.py): percent bias = aggregate ratio of means (ratio of
sums), plot-clustered bootstrap B = 2000, seed 20260725, percentile 2.5 / 97.5, plot
draws shared across arms within a cohort.

Positivity mask FIX applied: predicted side only (MCuFt_PRED > 0 and BA_PRED > 0).
The observed side is required present (notna) and non-negative, not strictly positive.

Cohorts:
  FULL  = corrected observed side, n = 3,053  (obs_cohort_volumes_corrected.csv)
  CLEAN = uncut, fully accounted subset, n = 2,525 (cohort_clean.csv + retain_looser.csv)

Configurations for the three-way comparison:
  C1 = bin-fixed_20260804  scorecards (finalrun_20260804/out_full),  FULL cohort
  C2 = bin-balfix_20260804 scorecards (balfix_run_20260804/out_full), FULL cohort
  C3 = bin-balfix_20260804 scorecards,                                CLEAN cohort
"""
import glob, os, sys, numpy as np, pandas as pd

BAL = "/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804/out_full"
FIX = "/fs/scratch/PUOM0008/crsfaaron/finalrun_20260804/out_full"
AUD = "/fs/scratch/PUOM0008/crsfaaron/obsaudit_20260805"
OUT = "/fs/scratch/PUOM0008/crsfaaron/balfix_gate_20260805"
OBC = AUD + "/obs_cohort_volumes_corrected.csv"
OB0 = "/fs/scratch/PUOM0008/crsfaaron/cull_test_20260728/obs_cohort_volumes.csv"
BFO = "/fs/scratch/PUOM0008/crsfaaron/bdft_gate_20260804/obs_cohort_bdft.csv"
SEED, B = 20260725, 2000
ARMS = ["fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate"]
GREG = ["organon", "conus_spdep", "conus_climate"]
NATIVE = ["fvs_base", "fvs_regional"]
os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 260)


def load_sc(d, tag):
    fs = sorted(glob.glob(d + "/scorecard_fin_armA_s*.csv"))
    sc = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    print("[%s] shards=%d rows=%d plots=%d models=%d"
          % (tag, len(fs), len(sc), sc.PLOT.nunique(), sc.model.nunique()))
    return sc, len(fs)


print("=" * 78)
print("PART 4a. SCORECARD SANITY GATES")
print("=" * 78)
sc, nsh = load_sc(BAL, "balfix")
scf, nshf = load_sc(FIX, "binfixed")

G = {}
G["rows_16605"] = (len(sc) == 16605)
G["plots_3321"] = (sc.PLOT.nunique() == 3321)
G["models_5"] = (sorted(sc.model.unique()) == sorted(ARMS))
ndup = int(sc.duplicated(["PLOT", "model"]).sum())
G["no_dup_plot_model"] = (ndup == 0)
G["shards_32"] = (nsh == 32)
print("GATE 1  rows==16605 : %s (%d)" % (G["rows_16605"], len(sc)))
print("GATE 1  plots==3321 : %s (%d)" % (G["plots_3321"], sc.PLOT.nunique()))
print("GATE 1  5 models    : %s (%s)" % (G["models_5"], sorted(sc.model.unique())))
print("GATE 1  dup (PLOT,model) pairs : %d  -> %s" % (ndup, G["no_dup_plot_model"]))
print("GATE 1  shards      : %d" % nsh)

# ---- GATE 4: fvs_base / fvs_regional bit-identity against bin-fixed, FULL scorecard
print("\n" + "-" * 78)
print("GATE 4. NATIVE-ARM BIT-IDENTITY, SCORECARD, ALL 3,321 PLOTS")
print("-" * 78)
num = ["BA_PRED", "TPA_PRED", "QMD_PRED", "MCuFt_PRED", "BdFt_PRED",
       "BA_OBS", "BA_OBS_WHOLE", "MCuFt_OBS", "BdFt_OBS", "N_COHORT",
       "PERIOD_YR", "N_SUM_ROWS", "YEAR_FIRST", "YEAR_LAST"]
num = [c for c in num if c in sc.columns and c in scf.columns]
ident = {}
for arm in ARMS:
    a = sc.loc[sc.model == arm].set_index("PLOT").sort_index()
    b = scf.loc[scf.model == arm].set_index("PLOT").sort_index()
    common = a.index.intersection(b.index)
    a = a.loc[common]; b = b.loc[common]
    mx = 0.0; worst = ""
    for c in num:
        dd = np.nanmax(np.abs(a[c].values.astype(float) - b[c].values.astype(float)))
        if not np.isfinite(dd):
            dd = 0.0
        if dd > mx:
            mx = dd; worst = c
    exact = all(a[c].astype(float).equals(b[c].astype(float)) for c in num)
    ident[arm] = mx
    print("  %-14s n=%d  maxabs diff = %.6g %s  exact_equal=%s"
          % (arm, len(common), mx, ("(" + worst + ")") if worst else "", exact))
G["native_bit_identical"] = all(ident[a] == 0.0 for a in NATIVE)
G["greg_moved"] = all(ident[a] > 0.0 for a in GREG)
print("GATE 4  fvs_base and fvs_regional bit-identical : %s" % G["native_bit_identical"])
print("        Greg arms moved (expected)              : %s" % G["greg_moved"])
if not G["native_bit_identical"]:
    print("\n*** STOP. A NATIVE ARM MOVED. LEAKAGE INTO THE NON-GREG PATH. ***")

# ---------------------------------------------------------------- cohorts
print("\n" + "=" * 78)
print("PART 2a. COHORT CONSTRUCTION AND THE POSITIVITY MASK FIX")
print("=" * 78)
obc = pd.read_csv(OBC)
ob0 = pd.read_csv(OB0)
bfo = pd.read_csv(BFO)
clean = set(pd.read_csv(AUD + "/cohort_clean.csv").PLOT.astype("int64"))
loose = set(pd.read_csv(AUD + "/retain_looser.csv").PLOT.astype("int64"))
CLEANSET = clean | loose
print("corrected observed file rows = %d   legacy = %d   clean+loose = %d"
      % (len(obc), len(ob0), len(CLEANSET)))
print("zero-observed gross plots in corrected file = %d"
      % int((obc.OBS_CFGRS == 0).sum()))


def build(scx, obx, keep=None, mask="new"):
    d = scx.merge(obx, on="PLOT", how="inner")
    if mask == "new":
        m = ((d.MCuFt_PRED > 0) & (d.BA_PRED > 0)
             & d.OBS_CFGRS.notna() & d.OBS_CFNET.notna() & d.OBS_BA.notna()
             & (d.OBS_CFGRS >= 0) & (d.OBS_CFNET >= 0) & (d.OBS_BA >= 0))
    else:
        m = ((d.MCuFt_PRED > 0) & (d.OBS_CFNET > 0) & (d.OBS_CFGRS > 0)
             & (d.BA_PRED > 0) & (d.OBS_BA > 0))
    d = d.loc[m].copy()
    if keep is not None:
        d = d.loc[d.PLOT.isin(keep)].copy()
    d["ARM"] = d["model"]
    cnt = d.groupby("PLOT").ARM.nunique()
    d = d.loc[d.PLOT.isin(set(cnt[cnt == 5].index))].copy()
    return d


# mask comparison, same scorecard and same observed file, old vs new mask
for lab, obx, keep in [("FULL", obc, None), ("CLEAN", obc, CLEANSET)]:
    dn = build(sc, obx, keep, "new")
    do = build(sc, obx, keep, "old")
    print("  %-6s cohort: NEW mask n=%d   OLD mask n=%d   delta=%+d"
          % (lab, dn.PLOT.nunique(), do.PLOT.nunique(),
             dn.PLOT.nunique() - do.PLOT.nunique()))

D_FULL = build(sc, obc, None, "new")
D_CLEAN = build(sc, obc, CLEANSET, "new")
D_FULL_FIX = build(scf, obc, None, "new")
D_CLEAN_FIX = build(scf, obc, CLEANSET, "new")
print("FULL  n = %d   CLEAN n = %d" % (D_FULL.PLOT.nunique(), D_CLEAN.PLOT.nunique()))
G["full_3053"] = (D_FULL.PLOT.nunique() == 3053)
G["clean_2525"] = (D_CLEAN.PLOT.nunique() == 2525)
print("cohort sizes match the audit (3053 / 2525): %s / %s"
      % (G["full_3053"], G["clean_2525"]))

# ---------------------------------------------------------------- gate engine
QTY = [("VOL_GROSS", "MCuFt_PRED", "OBS_CFGRS"),
       ("VOL_NET", "MCuFt_PRED", "OBS_CFNET"),
       ("VOL_SOUND", "MCuFt_PRED", "OBS_CFSND"),
       ("BA", "BA_PRED", "OBS_BA"),
       ("TPA", "TPA_PRED", "OBS_TPA")]


def gate(d, qtys, tag):
    plots = np.sort(d.PLOT.unique()); P = len(plots)
    o1 = d.drop_duplicates("PLOT").set_index("PLOT").loc[plots]
    obscols = sorted({oc for _, _, oc in qtys})
    predcols = sorted({pc for _, pc, _ in qtys})
    OBSM = np.column_stack([o1[c].values.astype(float) for c in obscols])
    PRED = {pc: np.zeros((P, len(ARMS))) for pc in predcols}
    for j, a in enumerate(ARMS):
        ga = d.loc[d.ARM == a].set_index("PLOT").loc[plots]
        for pc in predcols:
            PRED[pc][:, j] = ga[pc].values.astype(float)
    E2 = {}
    for lab, pc, oc in qtys:
        k = obscols.index(oc)
        E2[lab] = (PRED[pc] - OBSM[:, [k]]) ** 2
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, P, size=(B, P))
    W = np.empty((B, P), dtype=np.float64)
    for b in range(B):
        W[b] = np.bincount(idx[b], minlength=P)
    NUM = {pc: W @ PRED[pc] for pc in predcols}
    DO = W @ OBSM
    SQ = {lab: W @ E2[lab] for lab in E2}
    rows = []
    for lab, pc, oc in qtys:
        k = obscols.index(oc)
        for j, a in enumerate(ARMS):
            g = d.loc[d.ARM == a]
            p = g[pc].values.astype(float); o = g[oc].values.astype(float); e = p - o
            pt = 100.0 * (p.sum() - o.sum()) / o.sum()
            rep = 100.0 * (NUM[pc][:, j] / DO[:, k] - 1.0)
            lo, hi = np.percentile(rep, [2.5, 97.5])
            rr = 100.0 * np.sqrt(SQ[lab][:, j] / P) / (DO[:, k] / P)
            rlo, rhi = np.percentile(rr, [2.5, 97.5])
            rows.append(dict(CFG=tag, QTY=lab, ARM=a, n_plots=P, bias_pct=pt,
                             lo=lo, hi=hi, se=rep.std(ddof=1),
                             ratio=p.sum() / o.sum(),
                             ratio_lo=1.0 + lo / 100.0, ratio_hi=1.0 + hi / 100.0,
                             mean_obs=o.mean(), mean_pred=p.mean(),
                             MAE=float(np.mean(np.abs(e))),
                             RMSE=float(np.sqrt(np.mean(e ** 2))),
                             rRMSE_pct=float(100.0 * np.sqrt(np.mean(e ** 2)) / o.mean()),
                             rRMSE_lo=rlo, rRMSE_hi=rhi))
    return pd.DataFrame(rows)


print("\n" + "=" * 78)
print("PART 2b. CUBIC GATE")
print("=" * 78)
allr = []
for tag, dd in [("C1_binfixed_FULL", D_FULL_FIX),
                ("C2_balfix_FULL", D_FULL),
                ("C3_balfix_CLEAN", D_CLEAN),
                ("C4_binfixed_CLEAN", D_CLEAN_FIX)]:
    print("  running %s  n=%d" % (tag, dd.PLOT.nunique())); sys.stdout.flush()
    allr.append(gate(dd, QTY, tag))
cub = pd.concat(allr, ignore_index=True)
cub.to_csv(OUT + "/gate_balfix_cubic.csv", index=False)

# ---------------------------------------------------------------- board foot
print("\n" + "=" * 78)
print("PART 2c. BOARD FOOT GATE")
print("=" * 78)
BQ = [("BDFT_GROSS", "BdFt_PRED", "OBS_BFGRS"), ("BDFT_NET", "BdFt_PRED", "OBS_BFNET")]


def gate_bf(d, tag):
    d = d.merge(bfo[["PLOT", "OBS_BFGRS"]], on="PLOT", how="left")
    d = d.loc[(d.BdFt_PRED > 0) & (d.OBS_BFGRS > 0) & (d.OBS_BFNET > 0)].copy()
    cnt = d.groupby("PLOT").ARM.nunique()
    d = d.loc[d.PLOT.isin(set(cnt[cnt == 5].index))].copy()
    return gate(d, BQ, tag)


bfr = []
for tag, dd in [("C1_binfixed_FULL", D_FULL_FIX),
                ("C2_balfix_FULL", D_FULL),
                ("C3_balfix_CLEAN", D_CLEAN),
                ("C4_binfixed_CLEAN", D_CLEAN_FIX)]:
    r = gate_bf(dd, tag)
    print("  %s board-foot cohort n = %d" % (tag, int(r.n_plots.iloc[0]))); sys.stdout.flush()
    bfr.append(r)
bdf = pd.concat(bfr, ignore_index=True)
bdf.to_csv(OUT + "/gate_balfix_bdft.csv", index=False)

full = pd.concat([cub, bdf], ignore_index=True)
full.to_csv(OUT + "/gate_balfix_all.csv", index=False)

# ---------------------------------------------------------------- printout
def show(cfg):
    print("\n" + "=" * 78)
    print("CONFIG %s" % cfg)
    print("=" * 78)
    t = full.loc[full.CFG == cfg]
    print("  n plots (cubic) = %d   n plots (board foot) = %d"
          % (int(t.loc[t.QTY == "VOL_GROSS"].n_plots.iloc[0]),
             int(t.loc[t.QTY == "BDFT_GROSS"].n_plots.iloc[0])))
    hdr = ("%-15s %-11s %8s %18s %9s %19s" %
           ("ARM", "QTY", "bias%", "95% CI", "rRMSE%", "rRMSE 95% CI"))
    print(hdr); print("-" * len(hdr))
    for a in ARMS:
        for q in ["VOL_GROSS", "VOL_NET", "BDFT_GROSS", "BDFT_NET", "BA", "TPA"]:
            r = t.loc[(t.ARM == a) & (t.QTY == q)]
            if not len(r):
                continue
            r = r.iloc[0]
            print("%-15s %-11s %+8.2f  (%+7.2f, %+7.2f) %9.2f  (%7.2f, %7.2f)"
                  % (a, q, r.bias_pct, r.lo, r.hi, r.rRMSE_pct, r.rRMSE_lo, r.rRMSE_hi))
    print("\n  agreement ratios (predicted / observed, aggregate)")
    for a in ARMS:
        rb = t.loc[(t.ARM == a) & (t.QTY == "BA")].iloc[0]
        rt = t.loc[(t.ARM == a) & (t.QTY == "TPA")].iloc[0]
        print("    %-15s BA %.4f (%.4f, %.4f)   stem retention %.4f (%.4f, %.4f)"
              % (a, rb.ratio, rb.ratio_lo, rb.ratio_hi,
                 rt.ratio, rt.ratio_lo, rt.ratio_hi))


for c in ["C1_binfixed_FULL", "C2_balfix_FULL", "C3_balfix_CLEAN", "C4_binfixed_CLEAN"]:
    show(c)

print("\n" + "=" * 78)
print("PART 3. THREE-CONFIGURATION COMPARISON, GROSS CUBIC VOLUME BIAS (pct)")
print("=" * 78)
print("%-15s %26s %26s %26s" % ("ARM", "C1 binfixed FULL",
                                "C2 balfix FULL", "C3 balfix CLEAN"))
for a in ARMS:
    s = ""
    for c in ["C1_binfixed_FULL", "C2_balfix_FULL", "C3_balfix_CLEAN"]:
        r = full.loc[(full.CFG == c) & (full.QTY == "VOL_GROSS") & (full.ARM == a)].iloc[0]
        s += "  %+7.2f (%+6.2f,%+6.2f)" % (r.bias_pct, r.lo, r.hi)
    print("%-15s %s" % (a, s))


def ovl(a, c1, c2, q="VOL_GROSS"):
    x = full.loc[(full.CFG == c1) & (full.QTY == q) & (full.ARM == a)].iloc[0]
    y = full.loc[(full.CFG == c2) & (full.QTY == q) & (full.ARM == a)].iloc[0]
    return not (x.hi < y.lo or y.hi < x.lo)


print("\nnon-overlap tests, gross cubic")
for a in ARMS:
    print("  %-15s C1 vs C2 separated: %-5s    C2 vs C3 separated: %-5s"
          % (a, not ovl(a, "C1_binfixed_FULL", "C2_balfix_FULL"),
             not ovl(a, "C2_balfix_FULL", "C3_balfix_CLEAN")))

import json
with open(OUT + "/sanity_scorecard.json", "w") as f:
    json.dump({k: bool(v) for k, v in G.items()}, f, indent=1)
print("\nSANITY (scorecard side):", {k: bool(v) for k, v in G.items()})
print("\nWROTE:", OUT + "/gate_balfix_all.csv")
