#!/usr/bin/env python3
"""runoff_engine_runs.py  (2026-10-30)

Model run-off driver for the FVS engines (contestant A, fvs-modern; contestant C, unmodified USDA FVS)
and the input builder for the fvs-conus projector (contestant B), on the catalog stand sample.

Stand sample and stand construction are the catalog's, unchanged: catalog/conus_catalog_runs.py
sample_plots() on fia_variant_assignment.csv with --n-stands 300 --seed 20261028, plot CN recovered
from the public PLOT table by (STATECD, COUNTYCD, PLOT, INVYR1), condition trees only, null TPA_UNADJ
records dropped, TPA divided by CONDPROP_UNADJ, the same species relabel, traj id C%04d by sample order.
--first N keeps the first N stands of that sample (the pilot); the full run uses every stand.

--mode long   100 yr at a 5 yr step (the catalog protocol; NOAUTOES via FVS_NOAUTOES=1 for EM and IE).
              Output: <out>/runoff_long_<variant>_<arm>.csv, the catalog long schema without plt_cn.
--mode near   one cycle whose length is the plot's own remeasurement interval rounded to whole years
              (k_cycle from obs/nearterm_obs_frame.csv): NUMCYCLE 1 and TIMEINT 0 k in the keyfile.
              Output: <out>/runoff_near_<variant>_<arm>.csv, one row per stand and config.
--b-inputs DIR  also write the projector inputs for these exact stands: DIR/standinit_<VAR>.csv and
              DIR/treeinit_<var>/<ST>_FVS_TREEINIT_PLOT.csv (original FIA SPCD, DIA in, HT ft, CR pct,
              TREE_COUNT = TPA_UNADJ / CONDPROP_UNADJ). --configs none skips the engine runs.
The engine library comes from FVS_LIB_DIR (fvs-modern lib/ for A, runoff/upstream_build/lib for C).
No PLT_CN is written to any output. Public FIA plot coordinates only (the projector's covariate join).
"""
import os, sys, re, argparse, logging
import numpy as np, pandas as pd

CAT = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog"
sys.path.insert(0, CAT)
import conus_catalog_runs as ccr   # sample_plots, read_state_tables, REGION, FIPS2ST, unit constants

log = logging.getLogger("runoff")


def kf_patch(text, k):
    """One cycle of k years: replace the TIMEINT and NUMCYCLE records (10-column keyword fields)."""
    text, n1 = re.subn(r"(?m)^TIMEINT.*$", f"{'TIMEINT':<10}{0:>10d}{int(k):>10d}", text)
    text, n2 = re.subn(r"(?m)^NUMCYCLE.*$", f"{'NUMCYCLE':<10}{1:>10d}", text)
    if n1 != 1 or n2 != 1:
        raise RuntimeError(f"keyfile patch failed (TIMEINT {n1}, NUMCYCLE {n2})")
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("--engine-variant", default="", help="library, crosswalk and keyfile variant when it differs from the sample variant (C for ACD = upstream NE)")
    ap.add_argument("--arm", required=True, help="native or gompit (survival set by FVS_GOMPIT env)")
    ap.add_argument("--assign", default=os.path.join(CAT, "fia_variant_assignment.csv"))
    ap.add_argument("--fia-dir", default=os.environ.get("FIA_DATA_DIR", os.path.expanduser("~/fia_data")))
    ap.add_argument("--n-stands", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20261028)
    ap.add_argument("--first", type=int, default=0, help="keep the first N stands of the catalog sample (0 = all)")
    ap.add_argument("--traj", default="", help="comma list of catalog traj ids to keep (e.g. C0012,C0040)")
    ap.add_argument("--seed-tol", type=float, default=0.20)
    ap.add_argument("--mode", choices=["long", "near"], required=True)
    ap.add_argument("--obs", default="/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/obs/nearterm_obs_frame.csv")
    ap.add_argument("--configs", default="default,calibrated", help="comma list, or none")
    ap.add_argument("--b-inputs", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    configs = [] if a.configs in ("", "none") else a.configs.split(",")
    v = a.variant.lower()
    ev = (a.engine_variant or a.variant).lower()

    if configs:
        import bakuzis_uncertainty_comparison as buc
        from fia_stand_generator import build_stand_init, build_tree_init
        from inventory_keyfile import _variant_crosswalk
        import inventory_keyfile as _ik
        log.info("FVS_LIB_DIR=%s  lib=%s", buc.FVS_LIB_DIR, os.path.join(buc.FVS_LIB_DIR, f"FVS{ev}.so"))
        noautoes = os.environ.get("FVS_NOAUTOES", "") == "1"
        _orig_mk = _ik.make_inventory_keyfile
        _orig_tpl = buc.KEYFILE_TEMPLATE
        state = {"k": None}
        def _mk(*args, **kw):
            if noautoes:
                kw["calibration_keywords"] = (kw.get("calibration_keywords") or "") + "\nNOAUTOES"
            txt = _orig_mk(*args, **kw)
            return kf_patch(txt, state["k"]) if state["k"] else txt
        _ik.make_inventory_keyfile = _mk
        def set_cycle(k):
            # eastern (ne, acd) path formats buc.KEYFILE_TEMPLATE; the inventory path calls _mk above
            state["k"] = k
            buc.KEYFILE_TEMPLATE = kf_patch(_orig_tpl, k) if k else _orig_tpl
        xw = _variant_crosswalk(ev)
        def relabel(spcd):
            if spcd in xw: return spcd, False
            if spcd < 300 and 299 in xw: return 299, True
            if 998 in xw: return 998, True
            if 999 in xw: return 999, True
            return None, True
    else:
        noautoes = False

    assign = pd.read_csv(a.assign, dtype={"PLT_CN": str})
    smp = ccr.sample_plots(assign, v, a.n_stands, a.seed)
    if smp.empty:
        log.error("no fitting plots for %s", v); sys.exit(0)
    smp["traj"] = [f"C{i + 1:04d}" for i in range(len(smp))]     # catalog sid, by sample order
    if a.first > 0:
        smp = smp.iloc[: a.first].copy()
    if a.traj:
        smp = smp[smp["traj"].isin(a.traj.split(","))].copy()
    log.info("variant %s: %d stands (first=%d)", v, len(smp), a.first)

    obs = None
    if a.mode == "near":
        obs = pd.read_csv(a.obs).set_index("plot_key")

    cond, plot, tree = ccr.read_state_tables(a.fia_dir, sorted(smp["STATECD"].unique()))
    cond["PLT_CN"] = cond["PLT_CN"].astype(str); tree["PLT_CN"] = tree["PLT_CN"].astype(str)
    if not plot.empty: plot["CN"] = plot["CN"].astype(str)
    key = plot.set_index(["STATECD", "COUNTYCD", "PLOT", "INVYR"])["CN"]
    key = key[~key.index.duplicated(keep="first")]
    smp["PLT_CN"] = [key.get((int(r.STATECD), int(r.COUNTYCD), int(r.PLOT), int(r.INVYR1)), None) for r in smp.itertuples()]
    tree = tree[tree["STATUSCD"] == 1]

    rows, near, gate, b_si, b_ti = [], [], [], [], []
    for _, r in smp.iterrows():
        sid = r["traj"]; pk = r["plot_key"]
        base = dict(variant=v.upper(), traj=sid, mode=a.mode)
        if r["PLT_CN"] is None or pd.isna(r["PLT_CN"]):
            gate.append(dict(base, status="NO_PLOT_CN")); continue
        cn = str(r["PLT_CN"]); condid = int(r["CONDID"]) if pd.notna(r.get("CONDID", np.nan)) else 1
        c = cond[(cond["PLT_CN"] == cn) & (cond["CONDID"] == condid)]
        if c.empty: c = cond[cond["PLT_CN"] == cn]
        if c.empty:
            gate.append(dict(base, status="NO_COND")); continue
        c = c.iloc[0].copy()
        c = {kk: (None if (isinstance(val, float) and np.isnan(val)) else val) for kk, val in c.to_dict().items()}
        t = tree[(tree["PLT_CN"] == cn) & (tree["CONDID"] == condid)]
        t = t[t["TPA_UNADJ"].notna()]
        if t.empty:
            gate.append(dict(base, status="NO_TREES")); continue
        condprop = c.get("CONDPROP_UNADJ")
        condprop = float(condprop) if condprop is not None and float(condprop) > 0 else 1.0
        t = t.copy(); t["TPA_UNADJ"] = t["TPA_UNADJ"] / condprop
        balive = float(c.get("BALIVE") if c.get("BALIVE") is not None else np.nan)
        k = None
        if a.mode == "near":
            if pk not in obs.index or not bool(obs.loc[pk, "obs_ok"]):
                gate.append(dict(base, status="NO_OBS")); continue
            k = int(obs.loc[pk, "k_cycle"])
        lat = lon = np.nan
        if not plot.empty:
            p = plot[plot["CN"] == cn]
            if not p.empty: lat = float(p["LAT"].iloc[0]); lon = float(p["LON"].iloc[0])
        if a.b_inputs:
            st = int(c.get("STATECD", r["STATECD"]) or r["STATECD"]); cty = int(c.get("COUNTYCD", 0) or 0)
            b_si.append(dict(STAND_CN=sid, STAND_ID=f"{st:02d}00{cty:03d}_{sid}", STATE=st, INV_YEAR=int(r["INVYR1"]),
                             LATITUDE=lat, LONGITUDE=lon, plot_key=pk, K_CYCLE=k if k else np.nan, BALIVE=balive))
            b_ti.append(pd.DataFrame(dict(ST=ccr.FIPS2ST.get(st), STAND_CN=sid, TREE_COUNT=t["TPA_UNADJ"].values,
                                          SPECIES=t["SPCD"].astype(int).values, DIAMETER=t["DIA"].values,
                                          HT=t["HT"].values, CRRATIO=t["CR"].values)))
        if not configs:
            gate.append(dict(base, status="B_INPUT_ONLY", k_cycle=k)); continue
        sdf = build_stand_init(c, sid, ev)
        sdf["region"] = ccr.REGION.get(ev, 9)
        sdf["state"] = int(c.get("STATECD", r["STATECD"]) or r["STATECD"]); sdf["county"] = int(c.get("COUNTYCD", 0) or 0)
        if np.isfinite(lat): sdf["latitude"] = lat; sdf["longitude"] = lon
        sisp = c.get("SISP", np.nan)
        if sisp is not None and pd.notna(sisp) and int(sisp) > 0: sdf["site_species"] = int(sisp)
        t2 = t.copy(); rl = [relabel(int(s)) for s in t2["SPCD"]]
        t2["SPCD"] = [x[0] if x[0] is not None else -1 for x in rl]; t2 = t2[t2["SPCD"] > 0]
        if t2.empty:
            gate.append(dict(base, status="NO_MAPPED_TREES")); continue
        tdf = build_tree_init(t2, sid)
        set_cycle(k)
        for cfg in configs:
            try:
                s = buc.run_scenario(sdf, tdf, sid, ev, cfg)
            except Exception as e:
                gate.append(dict(base, config=cfg, status="RUN_FAIL", err=str(e)[:200])); continue
            if s is None or s.empty:
                gate.append(dict(base, config=cfg, status="EMPTY")); continue
            s = s.copy().sort_values("year")
            ba0 = float(s["atba"].iloc[0])
            seed_ok = bool(np.isfinite(balive) and balive > 0 and (abs(ba0 - balive) / balive <= a.seed_tol or abs(ba0 - balive) <= 1.0))
            arm = f"{cfg}_{a.arm}" + (f"({ev}_lib)" if ev != v else "")   # e.g. usfs_native(ne_lib) for C on ACD
            gate.append(dict(base, config=cfg, arm=arm, status="OK" if seed_ok else "SEED_FAIL", ba0_ft2ac=ba0, balive=balive,
                             tpa0=float(s["tpa"].iloc[0]), n_rows=len(s), k_cycle=k, noautoes=noautoes))
            if a.mode == "long":
                y0 = int(s["year"].iloc[0])
                rows.append(pd.DataFrame({"arm": arm, "group": v.upper(), "site": "all", "density": "all", "traj": sid,
                    "age": (s["year"] - y0).astype(float), "year": s["year"].astype(int),
                    "TPH": s["tpa"] * ccr.TPA_TPH, "BAPH": s["atba"] * ccr.FT2AC_M2HA,
                    "HT": s["attopht"] * ccr.FT_M, "VOL_NATIVE": s["tcuft"] * ccr.CUFTAC_M3HA,
                    "config": cfg, "seed_ok": seed_ok, "ba0_ft2ac": ba0, "balive": balive}))
            else:
                e = s.iloc[-1]; b = s.iloc[0]
                near.append(dict(variant=v.upper(), traj=sid, plot_key=pk, arm=arm, config=cfg, k_cycle=k,
                    years_run=int(e["year"] - b["year"]), seed_ok=seed_ok, balive=balive, ba0_ft2ac=ba0,
                    BA0=float(b["atba"]) * ccr.FT2AC_M2HA, TPH0=float(b["tpa"]) * ccr.TPA_TPH, HT0=float(b["attopht"]) * ccr.FT_M,
                    BA=float(e["atba"]) * ccr.FT2AC_M2HA, TPH=float(e["tpa"]) * ccr.TPA_TPH, HT=float(e["attopht"]) * ccr.FT_M,
                    VOL_NATIVE=float(e["tcuft"]) * ccr.CUFTAC_M3HA))
    os.makedirs(a.out, exist_ok=True)
    tag = f"{v}_{a.arm}"
    pd.DataFrame(gate).to_csv(os.path.join(a.out, f"runoff_{a.mode}_{tag}_gate.csv"), index=False)
    if rows:
        d = pd.concat(rows, ignore_index=True); d = d[d["TPH"] > 0].copy()
        d["QMD"] = np.sqrt(d["BAPH"] / d["TPH"] / (np.pi / 40000.0))
        d["VOL"] = 0.40 * d["BAPH"] * d["HT"]      # common volume, the scorer's form factor (fairness control)
        d.to_csv(os.path.join(a.out, f"runoff_long_{tag}.csv"), index=False)
    if near:
        d = pd.DataFrame(near)
        for p in ("", "0"):
            d["QMD" + p] = np.sqrt(d["BA" + p] / d["TPH" + p] / (np.pi / 40000.0))
        d["VOLC"] = 0.40 * d["BA"] * d["HT"]
        d.to_csv(os.path.join(a.out, f"runoff_near_{tag}.csv"), index=False)
    if a.b_inputs and b_si:
        os.makedirs(a.b_inputs, exist_ok=True)
        pd.DataFrame(b_si).to_csv(os.path.join(a.b_inputs, f"standinit_{v.upper()}.csv"), index=False)
        ti = pd.concat(b_ti, ignore_index=True); td = os.path.join(a.b_inputs, f"treeinit_{v}"); os.makedirs(td, exist_ok=True)
        for st, g in ti.groupby("ST"):
            g.drop(columns="ST").to_csv(os.path.join(td, f"{st}_FVS_TREEINIT_PLOT.csv"), index=False)
    g = pd.DataFrame(gate)
    log.info("status counts: %s", g["status"].value_counts().to_dict() if len(g) else {})
    print("RUNOFF_TASK_DONE", v, a.arm, a.mode, flush=True)


if __name__ == "__main__":
    main()
