#!/usr/bin/env python3
"""conus_catalog_runs.py  (2026-10-28)

Common-stand-sample engine runs for the CONUS configuration catalog.

For one FVS variant, draw up to --n-stands real FIA plots from the fitting data's own
variant assignment (fia_variant_assignment.csv, written from the CONUS remeasurement
pairs RDS on Cardinal), build FVS_StandInit / FVS_TreeInit rows from the per-state FIA
COND, PLOT and TREE tables, and project each stand 100 years at a 5 year step under the
fvs-modern engine with the default and the calibrated MAP parameter sets. The survival
model (native FVS mortality or GOMPIT) is set by the environment (FVS_GOMPIT=1 and
FVS_GOMPIT_COEF), exactly as in bakuzis_uncertainty_comparison.py, so one invocation is
one survival arm and the arm is written into the output rows from --arm.

The stand sample is the same for every parameter set and survival arm (seeded on the
variant), which is what makes the catalog arms comparable (plan red team, defect 1).
No posterior draws. No regeneration keywords. Public FIA plot coordinates only.

Output (long format read by bakuzis_stress_v2.R --source=long):
  <out>/catalog_<variant>_<arm>.csv
    arm, group (= variant), site ("all"), density ("all"), traj (= stand id), age,
    year, TPH, BAPH, QMD, HT, VOL, config, plt_cn, seed_ok, ba0_ft2ac, balive
  <out>/catalog_<variant>_<arm>_seedgate.csv  one row per stand and config
"""
import os, sys, argparse, logging
import numpy as np, pandas as pd

TPA_TPH = 2.4710538; FT2AC_M2HA = 0.2295684; FT_M = 0.3048; CUFTAC_M3HA = 0.0699727

# FVS variant -> USFS region code written to FVS_StandInit.region. [IMPLICIT: standard
# variant geography; FVS only uses region/forest for location-dependent defaults.]
REGION = {"acd": 9, "ne": 9, "ls": 9, "cs": 9, "sn": 8, "bm": 6, "ec": 6, "pn": 6, "wc": 6,
          "so": 6, "oc": 6, "op": 6, "ca": 5, "nc": 5, "ws": 5, "ie": 1, "em": 1, "kt": 1,
          "ci": 4, "tt": 4, "ut": 4, "cr": 2}

FIPS2ST = {1: "AL", 4: "AZ", 5: "AR", 6: "CA", 8: "CO", 9: "CT", 10: "DE", 12: "FL", 13: "GA",
           16: "ID", 17: "IL", 18: "IN", 19: "IA", 20: "KS", 21: "KY", 22: "LA", 23: "ME",
           24: "MD", 25: "MA", 26: "MI", 27: "MN", 28: "MS", 29: "MO", 30: "MT", 31: "NE",
           32: "NV", 33: "NH", 34: "NJ", 35: "NM", 36: "NY", 37: "NC", 38: "ND", 39: "OH",
           40: "OK", 41: "OR", 42: "PA", 44: "RI", 45: "SC", 46: "SD", 47: "TN", 48: "TX",
           49: "UT", 50: "VT", 51: "VA", 53: "WA", 54: "WV", 55: "WI", 56: "WY"}

log = logging.getLogger("catalog")


def sample_plots(assign, variant, n, seed):
    v = variant.upper()
    a = assign[assign["variant"] == v]
    if v == "ACD":   # Acadian: Northeast-variant plots in ME, NH, VT [ASSUMPTION]
        a = assign[(assign["variant"] == "NE") & assign["STATECD"].isin([23, 33, 50])]
    if len(a) == 0:
        return a
    if len(a) > n:
        a = a.sample(n=n, random_state=seed)
    return a.reset_index(drop=True)


def read_state_tables(fia_dir, states):
    cond, plot, tree = [], [], []
    for s in states:
        ab = FIPS2ST.get(int(s))
        if ab is None: continue
        pc = os.path.join(fia_dir, f"{ab}_COND.csv"); pp = os.path.join(fia_dir, f"{ab}_PLOT.csv"); pt = os.path.join(fia_dir, f"{ab}_TREE.csv")
        if not (os.path.exists(pc) and os.path.exists(pt)):
            log.warning("missing FIA tables for %s", ab); continue
        cond.append(pd.read_csv(pc, usecols=lambda c: c in ("PLT_CN", "CONDID", "COND_STATUS_CD", "SICOND", "SIBASE", "SISP", "BALIVE", "FORTYPCD", "STDAGE", "ASPECT", "SLOPE", "ELEV", "STATECD", "COUNTYCD", "CONDPROP_UNADJ"), low_memory=False))
        if os.path.exists(pp):
            plot.append(pd.read_csv(pp, usecols=lambda c: c in ("CN", "LAT", "LON", "STATECD", "COUNTYCD", "PLOT", "INVYR"), low_memory=False))
        tree.append(pd.read_csv(pt, usecols=lambda c: c in ("PLT_CN", "CONDID", "STATUSCD", "SPCD", "DIA", "HT", "CR", "TPA_UNADJ"), low_memory=False))
    return (pd.concat(cond, ignore_index=True) if cond else pd.DataFrame(),
            pd.concat(plot, ignore_index=True) if plot else pd.DataFrame(),
            pd.concat(tree, ignore_index=True) if tree else pd.DataFrame())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("--arm", required=True, help="native or gompit; the survival model is set by the environment")
    ap.add_argument("--assign", required=True, help="fia_variant_assignment.csv")
    ap.add_argument("--fia-dir", default=os.environ.get("FIA_DATA_DIR", os.path.expanduser("~/fia_data")))
    ap.add_argument("--n-stands", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20261028)
    ap.add_argument("--seed-tol", type=float, default=0.20)
    ap.add_argument("--out", required=True)
    ap.add_argument("--configs", default="default,calibrated")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    import bakuzis_uncertainty_comparison as buc        # run_scenario, paths, fvs2py setup
    from fia_stand_generator import build_stand_init, build_tree_init
    from inventory_keyfile import _variant_crosswalk

    v = a.variant.lower()
    # Species the compiled variant does not recognise are relabelled the way FVS itself does on input: softwoods
    # (SPCD < 300) to the variant's other-softwood code (FIA 299), hardwoods to other-hardwood (998), or to the
    # all-species code (999) where that is the only catch-all. Added October 29 after 165 of 300 SN stands, 76 BM
    # and 73 CS were lost to the strict crosswalk. Counts go into the seed gate file.
    xw = _variant_crosswalk(v)
    def relabel(spcd):
        if spcd in xw: return spcd, False
        if spcd < 300 and 299 in xw: return 299, True
        if 998 in xw: return 998, True
        if 999 in xw: return 999, True
        return None, True
    assign = pd.read_csv(a.assign, dtype={"PLT_CN": str})
    smp = sample_plots(assign, v, a.n_stands, a.seed)
    if smp.empty:
        log.error("no fitting plots assigned to variant %s; nothing to run", v); sys.exit(0)
    log.info("variant %s: %d stands sampled from %d assigned", v, len(smp), (assign["variant"] == v.upper()).sum())
    cond, plot, tree = read_state_tables(a.fia_dir, sorted(smp["STATECD"].unique()))
    cond["PLT_CN"] = cond["PLT_CN"].astype(str); tree["PLT_CN"] = tree["PLT_CN"].astype(str)
    if not plot.empty: plot["CN"] = plot["CN"].astype(str)
    # The pairs RDS carries PLT_CN as a double-cast (denormal) value, so the plot CN is recovered from the
    # public PLOT table by (STATECD, COUNTYCD, PLOT, INVYR), which the pairs file also carries.
    key = plot.set_index(["STATECD", "COUNTYCD", "PLOT", "INVYR"])["CN"]
    key = key[~key.index.duplicated(keep="first")]
    smp["PLT_CN"] = [key.get((int(r.STATECD), int(r.COUNTYCD), int(r.PLOT), int(r.INVYR1)), None) for r in smp.itertuples()]
    log.info("plot CN recovered for %d of %d sampled plots", smp["PLT_CN"].notna().sum(), len(smp))
    tree = tree[tree["STATUSCD"] == 1]

    rows, gate = [], []
    for i, r in smp.iterrows():
        if r["PLT_CN"] is None or pd.isna(r["PLT_CN"]):
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn="", status="NO_PLOT_CN")); continue
        cn = str(r["PLT_CN"]); condid = int(r["CONDID"]) if pd.notna(r.get("CONDID", np.nan)) else 1
        c = cond[(cond["PLT_CN"] == cn) & (cond["CONDID"] == condid)]
        if c.empty: c = cond[cond["PLT_CN"] == cn]
        if c.empty:
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, status="NO_COND")); continue
        c = c.iloc[0].copy()
        # FIA leaves STDAGE, ASPECT, SLOPE, ELEV and SICOND blank on some conditions; the generator's
        # defaults (40 yr, 0, 15 percent, 1000 ft, 60) apply, and the blank is recorded in the seed gate file.
        blanks = [k for k in ("STDAGE", "ASPECT", "SLOPE", "ELEV", "SICOND") if pd.isna(c.get(k, np.nan))]
        c = {k: (None if (isinstance(val, float) and np.isnan(val)) else val) for k, val in c.to_dict().items()}
        # Trees of this condition only (no all-plot fallback: that is what put SN and CS stands at 1.25 to 1.64 x
        # BALIVE on October 29), expanded to the condition acre by CONDPROP_UNADJ, which is how BALIVE is defined
        # (sum of BA x TPA_UNADJ / CONDPROP_UNADJ) and why partial conditions seeded at exactly 0.75 x before.
        t = tree[(tree["PLT_CN"] == cn) & (tree["CONDID"] == condid)]
        if t.empty:
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, status="NO_TREES")); continue
        condprop = c.get("CONDPROP_UNADJ")
        condprop = float(condprop) if condprop is not None and float(condprop) > 0 else 1.0
        t = t.copy(); t["TPA_UNADJ"] = t["TPA_UNADJ"] / condprop
        sid = f"C{i + 1:04d}"
        try:
            sdf = build_stand_init(c, sid, v)
        except Exception as e:
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, status="STAND_BUILD_FAIL", err=str(e)[:200])); continue
        sdf["region"] = REGION.get(v, 9)
        sdf["state"] = int(c.get("STATECD", r["STATECD"]) or r["STATECD"]); sdf["county"] = int(c.get("COUNTYCD", 0) or 0)
        if not plot.empty:
            p = plot[plot["CN"] == cn]
            if not p.empty:
                sdf["latitude"] = float(p["LAT"].iloc[0]); sdf["longitude"] = float(p["LON"].iloc[0])
        sisp = c.get("SISP", np.nan)
        if pd.notna(sisp) and int(sisp) > 0: sdf["site_species"] = int(sisp)
        t = t.copy()
        rl = [relabel(int(s)) for s in t["SPCD"]]
        t["SPCD"] = [r[0] if r[0] is not None else -1 for r in rl]
        n_relabel = int(sum(r[1] for r in rl)); n_drop = int((t["SPCD"] < 0).sum())
        ba_rel_share = float((t.loc[[r[1] for r in rl], "DIA"].pow(2) * t.loc[[r[1] for r in rl], "TPA_UNADJ"]).sum() / max((t["DIA"].pow(2) * t["TPA_UNADJ"]).sum(), 1e-9))
        t = t[t["SPCD"] > 0]
        if t.empty:
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, status="NO_MAPPED_TREES")); continue
        tdf = build_tree_init(t, sid)
        balive = float(c.get("BALIVE") if c.get("BALIVE") is not None else np.nan)
        for cfg in a.configs.split(","):
            try:
                s = buc.run_scenario(sdf, tdf, sid, v, cfg)
            except Exception as e:
                log.error("%s %s %s failed: %s", sid, cn, cfg, e); gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, config=cfg, status="RUN_FAIL", err=str(e)[:200])); continue
            if s is None or s.empty:
                gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, config=cfg, status="EMPTY")); continue
            s = s.copy().sort_values("year")
            ba0 = float(s["atba"].iloc[0]); tpa0 = float(s["tpa"].iloc[0])
            seed_ok = (np.isfinite(balive) and balive > 0 and abs(ba0 - balive) / balive <= a.seed_tol)
            gate.append(dict(variant=v.upper(), arm=a.arm, plt_cn=cn, config=cfg, status="OK" if seed_ok else "SEED_FAIL", ba0_ft2ac=ba0, balive=balive, tpa0=tpa0, stdage=c.get("STDAGE", np.nan), fortypcd=c.get("FORTYPCD", np.nan), sicond=c.get("SICOND", np.nan), blank_fields=";".join(blanks), n_relabel=n_relabel, n_drop=n_drop, ba_relabel_share=round(ba_rel_share, 4), condprop=condprop))
            y0 = int(s["year"].iloc[0])
            rows.append(pd.DataFrame({
                "arm": f"{cfg}_{a.arm}", "group": v.upper(), "site": "all", "density": "all", "traj": sid,
                "age": (s["year"] - y0).astype(float), "year": s["year"].astype(int),
                "TPH": s["tpa"] * TPA_TPH, "BAPH": s["atba"] * FT2AC_M2HA,
                "HT": s["attopht"] * FT_M if "attopht" in s else np.nan,
                "VOL": s["tcuft"] * CUFTAC_M3HA if "tcuft" in s else np.nan,
                "config": cfg, "plt_cn": cn, "seed_ok": seed_ok, "ba0_ft2ac": ba0, "balive": balive}))
        if (i + 1) % 25 == 0:
            log.info("  %d / %d stands done", i + 1, len(smp)); write_out(a, v, rows, gate, partial=True)
    write_out(a, v, rows, gate, partial=False)
    print("CATALOG_TASK_DONE", v, a.arm, flush=True)


def write_out(a, v, rows, gate, partial):
    """Checkpoint: the output files are rewritten every 25 stands so a wall-time kill still leaves the
    stands done so far (a 'partial' column says so; the final write clears it)."""
    os.makedirs(a.out, exist_ok=True)
    g = pd.DataFrame(gate); g["partial"] = partial
    tmp = os.path.join(a.out, f".catalog_{v}_{a.arm}_seedgate.tmp"); g.to_csv(tmp, index=False); os.replace(tmp, os.path.join(a.out, f"catalog_{v}_{a.arm}_seedgate.csv"))
    if rows:
        d = pd.concat(rows, ignore_index=True)
        d = d[d["TPH"] > 0].copy()
        d["QMD"] = np.sqrt(d["BAPH"] / d["TPH"] / (np.pi / 40000.0)); d["partial"] = partial
        tmp = os.path.join(a.out, f".catalog_{v}_{a.arm}.tmp"); d.to_csv(tmp, index=False); os.replace(tmp, os.path.join(a.out, f"catalog_{v}_{a.arm}.csv"))
        if not partial: log.info("wrote %d rows; seed gate OK %d of %d", len(d), int((g["status"] == "OK").sum()), len(g))
    elif not partial:
        log.error("no rows written for %s %s", v, a.arm)


if __name__ == "__main__":
    main()
