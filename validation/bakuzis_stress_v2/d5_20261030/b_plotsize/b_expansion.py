"""diagB: does per-plot-size expansion (SUBPPROP_UNADJ / MICRPROP_UNADJ / MACRPROP_UNADJ) explain the
multi-condition BALIVE residual after D3 (trees / CONDPROP_UNADJ)? D5 review rerun (2026-10-30): the
plot-size expansion is now the shipped fia_stand_generator.condition_expansion_divisor (1 percent plot-size
match, CONDPROP_UNADJ fallback), and the trees that fall back are counted per stand. Ensemble (NE, ACD, 360 stands, same
sampling as production) and catalog (4,748 stands, catalog sample_plots). FIA data server-side only;
no identifiers written: outputs carry variant/scenario/replicate or catalog traj only."""
import os, sys, warnings, logging
import numpy as np, pandas as pd
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
W = "/fs/scratch/PUOM0008/crsfaaron/wt-fix2-20261030"
B = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027"; OUT = B + "/fix_b_20261030/r1"
for p in (W, W + "/deployment/fvs2py", W + "/calibration/python", B + "/fix_b_20261030"):
    sys.path.insert(0, p)
import bakuzis_uncertainty_comparison as buc; logging.disable(logging.CRITICAL)
import fia_stand_generator as g
import conus_catalog_runs as ccr
from pathlib import Path
fia = Path(os.environ.get("FIA_DATA_DIR", str(Path.home() / "fia_data")))
CC = ("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "SUBPPROP_UNADJ", "MICRPROP_UNADJ", "MACRPROP_UNADJ", "BALIVE")

def tree_ba(t, c, how):
    t = t.copy(); ba = 0.005454154 * t["DIA"] ** 2 * t["TPA_UNADJ"]
    if how == "condprop":
        d = c["CONDPROP_UNADJ"]; return float((ba / d).sum()) if d and d > 0 else np.nan
    tpa = t["TPA_UNADJ"]
    if how == "threshold":   # the first D5 rule (above 50 micro, below 2 macro), kept to show nothing moved
        cls = np.where(tpa > 50, "MICRPROP_UNADJ", np.where(tpa < 2, "MACRPROP_UNADJ", "SUBPPROP_UNADJ"))
        den = np.array([c.get(k, np.nan) for k in cls], dtype=float)
        den = np.where(np.isfinite(den) & (den > 0), den, c["CONDPROP_UNADJ"])
        return float((ba / den).sum())
    den, n_no_size, n_no_prop = g.condition_expansion_divisor(tpa.to_numpy(), c, return_counts=True)
    return float((ba / den).sum()), n_no_size, n_no_prop

def evaluate(stands, trees, cond_all, keycols):
    rows = []
    for key, cn, cid, balive in stands:
        cc = cond_all[(cond_all["PLT_CN"] == cn) & (cond_all["CONDID"] == cid)]
        if cc.empty: continue
        c = cc.iloc[0].to_dict()
        t = trees[(trees["PLT_CN"] == cn) & (trees["CONDID"] == cid)]
        ncond = int((cond_all["PLT_CN"] == cn).sum())
        r = dict(zip(keycols, key)); r.update(n_cond=ncond, balive=balive, condprop=c["CONDPROP_UNADJ"], subpprop=c.get("SUBPPROP_UNADJ"), micrprop=c.get("MICRPROP_UNADJ"),
            macrprop=c.get("MACRPROP_UNADJ"), n_micro=int((t["TPA_UNADJ"] > 50).sum()), n_macro=int((t["TPA_UNADJ"] < 2).sum()),
            ba_condprop=tree_ba(t, c, "condprop"), ba_threshold=tree_ba(t, c, "threshold"))
        r["ba_plotsize"], r["n_no_plot_size"], r["n_no_proportion"] = tree_ba(t, c, "plotsize")
        rows.append(r)
    d = pd.DataFrame(rows)
    d["r_condprop"] = d["ba_condprop"] / d["balive"]; d["r_plotsize"] = d["ba_plotsize"] / d["balive"]
    d["ba_change_pct"] = 100 * (d["ba_plotsize"] / d["ba_condprop"] - 1)
    d["shipped_vs_threshold_pct"] = 100 * (d["ba_plotsize"] / d["ba_threshold"] - 1)
    return d

def load(states, cns):
    cond, tfr = [], []
    for ab in states:
        cc = g._state_csv(ab, "COND", fia); tc = g._state_csv(ab, "TREE", fia)
        if cc is None or tc is None: continue
        c = pd.read_csv(cc, usecols=lambda x: x in CC, low_memory=False); c = c[c["PLT_CN"].isin(cns)]; cond.append(c)
        for ch in pd.read_csv(tc, usecols=lambda x: x in ("PLT_CN", "CONDID", "STATUSCD", "DIA", "TPA_UNADJ"), chunksize=2_000_000, low_memory=False):
            ch = ch[ch["PLT_CN"].isin(cns) & (ch["STATUSCD"] == 1) & ch["TPA_UNADJ"].notna()]
            if len(ch): tfr.append(ch)
    return pd.concat(cond, ignore_index=True), pd.concat(tfr, ignore_index=True)

def summarize(d, label):
    m = d[d["n_cond"] > 1]; dev = lambda r: (r - 1).abs() > 0.02
    print(f"\n== {label}: stands {len(d)}, multi-condition {len(m)}")
    for col in ("r_condprop", "r_plotsize"):
        print(f"  {col}: multi-cond median {m[col].median():.4f} q05 {m[col].quantile(.05):.4f} q95 {m[col].quantile(.95):.4f} min {m[col].min():.3f} max {m[col].max():.3f}; "
              f"off by >2%: {int(dev(m[col]).sum())}; single-cond off by >2%: {int(dev(d.loc[d['n_cond'] == 1, col]).sum())}")
    print(f"  stands whose seeded BA moves >2% when switching to plot-size expansion: {int((d['ba_change_pct'].abs() > 2).sum())} of {len(d)} ({100 * (d['ba_change_pct'].abs() > 2).mean():.2f}%); max |change| {d['ba_change_pct'].abs().max():.2f}%")
    print(f"  fallback trees: no plot size {int(d['n_no_plot_size'].sum())} in {int((d['n_no_plot_size'] > 0).sum())} stands "
          f"({int((m['n_no_plot_size'] > 0).sum())} multi-cond); no proportion {int(d['n_no_proportion'].sum())} in {int((d['n_no_proportion'] > 0).sum())} stands "
          f"({int((m['n_no_proportion'] > 0).sum())} multi-cond); stands whose BA differs shipped vs threshold rule by >1e-7 percent: "
          f"{int((d['shipped_vs_threshold_pct'].abs() > 1e-7).sum())} (multi-cond {int((m['shipped_vs_threshold_pct'].abs() > 1e-7).sum())}), max |diff| {d['shipped_vs_threshold_pct'].abs().max():.3g}%")
    print(f"  multi-cond r_plotsize exactly 1 within 1e-4: {int(((m['r_plotsize'] - 1).abs() < 1e-4).sum())} of {len(m)}")
    print(f"  multi-cond stands where SUBPPROP != CONDPROP: {int((m['subpprop'] != m['condprop']).sum())}; MICRPROP != CONDPROP: {int((m['micrprop'] != m['condprop']).sum())}")

# ---- ensemble stands (NE, ACD), production sampling
rows = []
for v in ("ne", "acd"):
    frames = []
    states = [g._state_abbrev(c) for c in g.VARIANT_STATES.get(v, ())]
    for ab in states:
        cc = g._state_csv(ab, "COND", fia)
        if cc is None: continue
        df = pd.read_csv(cc, usecols=lambda c: c in ("STATECD", "COUNTYCD", "PLOT", "INVYR", "COND_STATUS_CD", "SICOND", "SIBASE", "BALIVE", "FORTYPCD", "STDAGE", "ASPECT", "SLOPE", "ELEV", "PLT_CN", "CONDID", "CONDPROP_UNADJ"), low_memory=False)
        frames.append(df.loc[df["COND_STATUS_CD"] == 1].copy())
    cond = pd.concat(frames, ignore_index=True)
    sg = buc.get_species_groups(v); scen = buc.enumerate_scenarios(sg); picked = []
    for (sp, site, dens, num) in scen:
        s = g.stratified_sample_conditions(cond, site, dens, n_plots=5, seed=42 + num)
        for rep, (_, r) in enumerate(s.iterrows(), start=1):
            picked.append(((v.upper(), num, rep), int(r["PLT_CN"]), int(r["CONDID"]), float(r["BALIVE"])))
    call, tr = load(states, {p[1] for p in picked})
    rows.append(evaluate(picked, tr, call, ("variant", "scenario", "replicate")))
E = pd.concat(rows, ignore_index=True); E.drop(columns=[]).to_csv(OUT + "/ensemble_expansion_check.csv", index=False)
summarize(E, "ensemble NE+ACD")
for v in ("NE", "ACD"): summarize(E[E["variant"] == v], f"ensemble {v}")

# ---- catalog sample (4,748 stands)
asg = pd.read_csv(B + "/catalog/fia_variant_assignment.csv", dtype={"PLT_CN": str})
crow = []
for v in "acd bm ca ci cr cs ec em ie ls ne pn sn so tt ut wc ws".split():
    smp = ccr.sample_plots(asg, v, 300, 20261028); smp["traj"] = [f"C{i + 1:04d}" for i in range(len(smp))]
    states = [ccr.FIPS2ST[int(s)] for s in sorted(smp["STATECD"].unique()) if int(s) in ccr.FIPS2ST]
    plots = []
    for ab in states:
        pp = g._state_csv(ab, "PLOT", fia)
        if pp is None: continue
        plots.append(pd.read_csv(pp, usecols=lambda x: x in ("CN", "STATECD", "COUNTYCD", "PLOT", "INVYR"), low_memory=False))
    pl = pd.concat(plots, ignore_index=True); key = pl.drop_duplicates(["STATECD", "COUNTYCD", "PLOT", "INVYR"]).set_index(["STATECD", "COUNTYCD", "PLOT", "INVYR"])["CN"]
    smp["cn"] = [key.get((int(r.STATECD), int(r.COUNTYCD), int(r.PLOT), int(r.INVYR1)), np.nan) for r in smp.itertuples()]
    smp = smp[smp["cn"].notna()]
    call, tr = load(states, set(smp["cn"].astype(np.int64)))
    bal = call.set_index(["PLT_CN", "CONDID"])["BALIVE"]
    st = [((v.upper(), r.traj), int(r.cn), int(r.CONDID) if pd.notna(r.CONDID) else 1, float(bal.get((int(r.cn), int(r.CONDID) if pd.notna(r.CONDID) else 1), np.nan))) for r in smp.itertuples()]
    crow.append(evaluate(st, tr, call, ("variant", "traj")))
    print("catalog", v, "done", flush=True)
C = pd.concat(crow, ignore_index=True); C.to_csv(OUT + "/catalog_expansion_check.csv", index=False)
summarize(C, "catalog all variants")
print(C.assign(moved=C["ba_change_pct"].abs() > 2).groupby("variant").agg(n=("moved", "size"), multi=("n_cond", lambda x: int((x > 1).sum())), moved=("moved", "sum"), max_abs_change=("ba_change_pct", lambda x: round(x.abs().max(), 2))).to_string())
print("DIAGB_DONE")
