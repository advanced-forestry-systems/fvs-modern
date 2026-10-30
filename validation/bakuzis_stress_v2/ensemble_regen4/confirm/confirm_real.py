"""Real-stand confirmations for the 2026-10-30 fixes (D1, D2, D3). No plot identifiers are written.

D1  FVS keyword echo of one NE and one ACD posterior-draw run (draw 481, the first of the
    ensemble's 100 draw indices): per-species MAX SDI and multipliers as FVS read them,
    against the values the engine intended.
D2  MAP arm with the fixed overlay (NA species skipped) vs the as-run overlay (NA -> 0) vs
    keywords only, on three stands per variant; and the draw arm with no overlay.
D3  BALIVE reproduction for every ensemble stand (360), by multi-condition status.
"""
import os, re, sys, glob, json, sqlite3, warnings, logging
import numpy as np, pandas as pd
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
F = "/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030"
X = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_20261030"
AUD = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/audit_20261030"
for p in (F, F + "/deployment/fvs2py", F + "/calibration/python"):
    sys.path.insert(0, p)
import bakuzis_uncertainty_comparison as buc
logging.disable(logging.CRITICAL)
import fia_stand_generator as g
from config.uncertainty import UncertaintyEngine
from config.config_loader import FvsConfigLoader
import config.config_loader as cl
from fvs2py import FVS
from pathlib import Path
CONFIG_DIR = F + "/config"; LIB = F + "/lib"
OUT = X + "/confirm"; os.makedirs(OUT, exist_ok=True)
fia_dir = Path(os.environ.get("FIA_DATA_DIR", str(Path.home() / "fia_data")))
which = sys.argv[1] if len(sys.argv) > 1 else "all"

# ---------------------------------------------------------------- D3
if which in ("all", "d3"):
    d3 = pd.read_csv(AUD + "/d3_stands.csv")[["variant", "scenario", "replicate", "multi", "contaminated", "balive", "ens_ba0"]]
    rows = []
    for v in ("ne", "acd"):
        cond = g.load_fia_conditions(v, fia_dir)
        scen = buc.enumerate_scenarios(buc.get_species_groups(v))
        picked = []
        for (sp, site, dens, num) in scen:
            s = g.stratified_sample_conditions(cond, site, dens, n_plots=5, seed=42 + num)
            for rep, (_, r) in enumerate(s.iterrows(), start=1):
                picked.append((num, rep, r))
        trees = g.load_fia_trees(v, sorted({int(p[2]["PLT_CN"]) for p in picked}), fia_dir)
        for num, rep, r in picked:
            t = g.condition_trees(trees, r)
            ti = g.build_tree_init(t, "X")
            ba = float((0.005454154 * ti["diameter"] ** 2 * ti["tree_count"]).sum()) if len(ti) else 0.0
            rows.append(dict(variant=v.upper(), scenario=num, replicate=rep, balive_cond=float(r["BALIVE"]),
                             condprop=float(r.get("CONDPROP_UNADJ", np.nan)), ba_seeded_new=ba))
    d = pd.DataFrame(rows).merge(d3, on=["variant", "scenario", "replicate"], how="left")
    assert np.allclose(d["balive_cond"], d["balive"]), "sampled conditions differ from the audited ensemble stands"
    d["ratio_new"] = d["ba_seeded_new"] / d["balive"]
    d["ratio_old_ens"] = d["ens_ba0"] / d["balive"]
    d.drop(columns=["balive_cond"]).to_csv(OUT + "/d3_balive_reproduction.csv", index=False)
    print("D3 stands", len(d), "multi-condition", int(d["multi"].sum()))
    for lab, sub in (("multi", d[d["multi"]]), ("single", d[~d["multi"]]), ("all", d)):
        q = sub["ratio_new"].quantile([.05, .5, .95]).round(4).tolist()
        qo = sub["ratio_old_ens"].quantile([.05, .5, .95]).round(4).tolist()
        print(f"D3 {lab:6s} n={len(sub)} new seed/BALIVE q05,med,q95={q} within1pct={int((abs(sub['ratio_new'] - 1) <= .01).sum())} | old ens_full3 year0/BALIVE q05,med,q95={qo}")

# ---------------------------------------------------------------- D1 / D2
def old_apply_to_fvs(self, fvs_instance):
    """The as-run (pre-2026-10-30) SDI overlay: 'NA' -> 0, pad with 0."""
    sdi = self._find_sdi_param(self.config.get("categories", {}))
    if sdi is None:
        return {}
    arr = np.array([0.0 if isinstance(v, str) or v is None else v for v in sdi], dtype=np.float64)
    maxsp = fvs_instance.dims.get("maxspecies", len(arr))
    arr = np.pad(arr, (0, max(0, maxsp - len(arr))))[:maxsp]
    fvs_instance.set_species_attr("spsdi", arr)
    return {"spsdi": True}

NEW_APPLY = FvsConfigLoader.apply_to_fvs


def run_case(v, sdf, tdf, sid, blk, version, tag, overlay="new"):
    d = os.path.join(OUT, "runs", v, tag); os.makedirs(d, exist_ok=True)
    for f in glob.glob(d + "/*"): os.remove(f)
    db = os.path.join(d, "FVS_Data.db")
    c = sqlite3.connect(db); sdf.to_sql("fvs_standinit", c, index=False); tdf.to_sql("fvs_treeinit", c, index=False); c.close()
    kp = os.path.join(d, f"{v}_{sid}.key")
    open(kp, "w").write(buc.KEYFILE_TEMPLATE.format(stand_id=sid, db_path=db, calibration_keywords=blk or "** DEFAULT PARAMETERS"))
    FvsConfigLoader.apply_to_fvs = old_apply_to_fvs if overlay == "old" else NEW_APPLY
    try:
        fvs = FVS(lib_path=f"{LIB}/FVS{v}.so", config_version=version, config_dir=CONFIG_DIR)
        fvs.load_keyfile(kp); fvs.run()
        s = fvs.summary.copy().sort_values("year")
        sp = np.asarray(fvs.get_species_attr("spsdi"), dtype=float).copy()
    finally:
        FvsConfigLoader.apply_to_fvs = NEW_APPLY
    echo = open(glob.glob(d + "/*.out")[0], errors="ignore").read()
    return s, sp, echo


def echo_table(echo):
    sdi = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"\(CODE=\s*(\d+)\); MAX SDI=\s*([\d.]+)", echo)}
    mult = {}
    for m in re.finditer(r"^(BAIMULT|MORTMULT|HTGMULT)\s+DATE/CYCLE=\s*\d+; SPECIES=\s*\S*\s*\(CODE=\s*(\d+)\); MULTIPLIER=\s*([\d.]+)", echo, re.M):
        mult[(m.group(1), int(m.group(2)))] = float(m.group(3))
    errs = len(re.findall(r"FVS0\d ERROR|KEYWORD IGNORED", echo))
    return sdi, mult, errs


if which in ("all", "d12"):
    run_rows, echo_rows = [], []
    for v in ("ne", "acd"):
        eng = UncertaintyEngine(v, config_dir=CONFIG_DIR, seed=42)
        idx = list(np.random.default_rng(42).choice(eng.n_draws, size=100, replace=False))
        di = int(idx[0])
        draw = eng.get_draw(di)
        mult = eng.compute_multipliers_for_draw(draw)
        dkw = eng.generate_keywords_for_draw(draw, draw_idx=di, include_comments=False)
        open(f"{OUT}/{v}_draw{di}_fixed.kw", "w").write(dkw)
        cal = FvsConfigLoader(v, version="calibrated", config_dir=CONFIG_DIR)
        ckw = cal.generate_keywords(include_comments=False)
        sdi_map = cal._find_sdi_param(cal.config["categories"])
        scen = buc.enumerate_scenarios(buc.get_species_groups(v))
        for k, (sp, site, dens, num) in enumerate([scen[i] for i in (1, 13, 31)]):
            st = g.generate_real_stand(variant=v, site_class=site, density_class=dens, n_plots=5, seed=42 + num)
            if not st: continue
            sdf, tdf, _ = st[0]; sid = f"BK{num:02d}R1"
            sdf = sdf.copy(); tdf = tdf.copy(); sdf["stand_id"] = sid; tdf["stand_id"] = sid
            cases = [("default", None, "default", "new"), ("calmap_asrun_overlay", ckw, "calibrated", "old"),
                     ("calmap_fixed_overlay", ckw, "calibrated", "new"), ("calmap_kwonly", ckw, "default", "new")]
            if k == 0:
                cases += [(f"draw{di}_fixed", dkw, "default", "new")]
            for tag, blk, ver, ov in cases:
                s, spsdi, echo = run_case(v, sdf, tdf, sid, blk, ver, f"s{num}_{tag}", ov)
                r0, r1 = s.iloc[0], s.iloc[-1]
                se, me, errs = echo_table(echo)
                run_rows.append(dict(variant=v.upper(), scenario=num, case=tag, y0=int(r0.year), yN=int(r1.year), ba0=r0.atba, baN=r1.atba,
                                     tpaN=r1.tpa, tcuftN=r1.tcuft, spsdi_zero=int((spsdi == 0).sum()), spsdi_min=float(spsdi.min()),
                                     n_sdimax_echo=len(se), n_mult_echo=len(me), fvs_errors=errs))
                if tag.startswith("draw"):
                    for i in range(len(spsdi)):
                        want = mult.get("sdi_values", np.full(len(spsdi), np.nan))
                        w = float(want[i]) if i < len(want) and np.isfinite(want[i]) else np.nan
                        echo_rows.append(dict(variant=v.upper(), fvs_sp=i + 1, kind="SDIMAX", intended=(round(w, 1) if np.isfinite(w) else np.nan),
                                              echoed=se.get(i + 1, np.nan), map=(np.nan if sdi_map is None or i >= len(sdi_map) or isinstance(sdi_map[i], str) else sdi_map[i]),
                                              fvs_spsdi_after_run=spsdi[i]))
                    for (kw, spn), val in me.items():
                        key = {"BAIMULT": "growth_mult", "MORTMULT": "mort_mult", "HTGMULT": "hg_mult"}[kw]
                        echo_rows.append(dict(variant=v.upper(), fvs_sp=spn, kind=kw, intended=round(float(mult[key][spn - 1]), 4), echoed=val))
    rr = pd.DataFrame(run_rows); rr.to_csv(OUT + "/d12_runs.csv", index=False)
    er = pd.DataFrame(echo_rows); er.to_csv(OUT + "/d1_echo_decode.csv", index=False)
    pd.set_option("display.width", 220)
    print(rr.round(2).to_string())
    for (v, kind), sub in er.groupby(["variant", "kind"]):
        e = sub.dropna(subset=["echoed"])
        tol = {"SDIMAX": 0.6, "MORTMULT": 0.0051}.get(kind, 1e-4 + 5e-5)
        ok = (abs(e["echoed"] - e["intended"]) <= tol).sum()
        print(f"D1 echo {v} {kind}: intended records={int(sub['intended'].notna().sum())} echoed={len(e)} match={ok} "
              f"distinct echoed values={e['echoed'].nunique()} range=({e['echoed'].min()}, {e['echoed'].max()})")
print("CONFIRM_DONE")
