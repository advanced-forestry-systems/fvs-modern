#!/usr/bin/env python3
"""
bakuzis_run.py -- long-horizon stratified FVS NE projection for the modified
Bakuzis matrix (Leary 1997) biological-realism assessment.

Adapted from run_driver_final.py (imported as a module, NOT modified).
Differences from the benchmark driver:
  * multi-cycle: TIMEINT 0 <CYCLE_LEN> / NUMCYCLE <N_CYC>  (default 10 yr x 15)
  * no remeasurement T2 side at all: this is a pure forward projection
  * plots stratified on SICOND (site) x initial TPA (density)
  * full FVS_Summary2 trajectory retained, not just first/last row
  * treelist output suppressed (15 cycles x 440 plots x 3 arms is too much I/O)

Usage:
  python3 bakuzis_run.py --stage select
  python3 bakuzis_run.py --stage run --models fvs_base,fvs_regional,organon \
        --nproc 32 --outdir out
"""
from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_driver_final as R  # noqa: E402

# --- BALFIX RERUN 2026-08-05: frozen PTBALT-repaired binary -----------------
BALFIX_BIN = "/users/PUOM0008/crsfaaron/fvs-modern/bin-balfix_20260804/FVSne"
BALFIX_MD5 = "b85e4c06995382425b01c756811d6d34"
import hashlib
_h = hashlib.md5(open(BALFIX_BIN, "rb").read()).hexdigest()
if _h != BALFIX_MD5:
    raise SystemExit("BINARY MD5 MISMATCH: got %s" % _h)
R.VARIANT_BINS["ne"] = BALFIX_BIN
print("[balfix] binary pinned: %s md5=%s" % (BALFIX_BIN, _h))

VARIANT = "ne"
CYCLE_LEN = 10
N_CYC = 15          # 150 years total
START_YEAR = 2000   # nominal; FVS Age comes from the standinit age field

# --- monkeypatch the cycle keyword for a long multi-cycle projection --------
def _long_cycle_kw(ncyc: int) -> str:
    return ("TIMEINT            0{c:>10d}\n"
            "NUMCYCLE  {n:>10d}").format(c=CYCLE_LEN, n=N_CYC)


R.build_cycle_kw = _long_cycle_kw

# --- monkeypatch the keyfile: drop the treelist, keep the summary ----------
R.KEYFILE = """\
STDIDENT
{sid}
DATABASE
DSNIN
{db}
DSNOUT
{db}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'
ENDSQL
END
DATABASE
SUMMARY            2
END
FIAVBC
{cycle_kw}
{htg_kw}
{calib_kw}
PROCESS
STOP
"""

# ---------------------------------------------------------------------------
# Stage 1: stratified plot selection
# ---------------------------------------------------------------------------
AGE_LO, AGE_HI = 20, 70     # initial stand age band, keeps a common age window
N_PER_CELL = 55             # 4 site classes x 2 density classes x 55 = 440


def select_plots(plot_list_csv: str, outcsv: str, seed: int = 20260804) -> None:
    states = R.VARIANTS[VARIANT]["states"]
    keep = set(pd.read_csv(plot_list_csv, usecols=["PLOT"])["PLOT"]
               .astype("int64").tolist())
    print(f"[select] benchmark T1 plot pool: {len(keep)}")

    conds = R.load_cond_for_cns(states, keep)
    trees = R.load_trees_for_cns(states, keep)
    print(f"[select] loaded {len(conds)} cond recs, {len(trees)} tree recs")

    # per-plot initial cohort stats (DIA >= 5 in, the FVS input list)
    trees["DIA"] = pd.to_numeric(trees["DIA"], errors="coerce")
    trees["TPA_UNADJ"] = pd.to_numeric(trees["TPA_UNADJ"], errors="coerce")
    t1 = trees.loc[trees["DIA"].notna() & (trees["DIA"] >= 5.0)
                   & (trees["STATUSCD"] == 1)].copy()
    t1["BA_TREE"] = 0.005454154 * t1["DIA"] ** 2 * t1["TPA_UNADJ"]
    agg = (t1.groupby("PLT_CN")
             .agg(TPA0=("TPA_UNADJ", "sum"), BA0=("BA_TREE", "sum"),
                  NTREE=("DIA", "size"))
             .reset_index())
    agg["PLT_CN"] = agg["PLT_CN"].astype("int64")

    c = conds.copy()
    c["PLT_CN"] = c["PLT_CN"].astype("int64")
    c = c.drop_duplicates(subset=["PLT_CN"], keep="first")
    want = [x for x in ["PLT_CN", "SICOND", "STDAGE", "STATECD", "SLOPE",
                        "ELEV", "FORTYPCD"] if x in c.columns]
    print(f"[select] cond columns available: {want}")
    d = agg.merge(c[want], on="PLT_CN", how="inner")
    for miss in ("STATECD",):
        if miss not in d.columns:
            d[miss] = 0

    d["SICOND"] = pd.to_numeric(d["SICOND"], errors="coerce")
    d["STDAGE"] = pd.to_numeric(d["STDAGE"], errors="coerce")
    n0 = len(d)
    d = d.loc[d["SICOND"].notna() & (d["SICOND"] >= 20) & (d["SICOND"] <= 110)]
    d = d.loc[d["STDAGE"].notna() & d["STDAGE"].between(AGE_LO, AGE_HI)]
    d = d.loc[d["NTREE"] >= 8]
    print(f"[select] after SICOND/STDAGE/NTREE filters: {n0} -> {len(d)}")

    # site class = SICOND quartiles (1 = poorest .. 4 = best)
    d["site_cls"] = pd.qcut(d["SICOND"], 4, labels=[1, 2, 3, 4]).astype(int)
    # density class = median split of initial TPA within site class
    d["dens_cls"] = (d.groupby("site_cls")["TPA0"]
                      .transform(lambda s: pd.qcut(s, 2, labels=["low", "high"]))
                      .astype(str))

    rng = np.random.default_rng(seed)
    parts = []
    for (s, dc), g in d.groupby(["site_cls", "dens_cls"]):
        n = min(N_PER_CELL, len(g))
        idx = rng.choice(len(g), size=n, replace=False)
        parts.append(g.iloc[idx])
        print(f"[select]   site {s} dens {dc}: pool {len(g)} -> take {n} "
              f"(SI {g['SICOND'].mean():.1f}, TPA0 {g['TPA0'].mean():.0f}, "
              f"age0 {g['STDAGE'].mean():.0f})")
    sel = pd.concat(parts, ignore_index=True)
    sel.to_csv(outcsv, index=False)
    print(f"[select] wrote {len(sel)} plots -> {outcsv}")
    print(sel.groupby("site_cls")["SICOND"].describe()[["min", "mean", "max"]])


# ---------------------------------------------------------------------------
# Stage 2: projection
# ---------------------------------------------------------------------------
_G = {}


def _load_shared(sel_path, states):
    """Load ONCE in the parent. Workers inherit via fork copy-on-write.

    Loading inside a Pool initializer made all 32 workers read the full
    per-state FIA TREE CSVs concurrently, which OOM-killed the job.
    """
    keep = set(pd.read_csv(sel_path)["PLT_CN"].astype("int64").tolist())
    _G["sel"] = pd.read_csv(sel_path).set_index("PLT_CN")
    _G["trees"] = R.load_trees_for_cns(states, keep)
    _G["conds"] = R.load_cond_for_cns(states, keep)
    _G["ecomap"] = R.load_ecoregion_map(keep)
    print(f"[run] shared load: {len(_G['trees'])} tree recs, "
          f"{len(_G['conds'])} cond recs, {len(_G['ecomap'])} ecoregions")
    sys.stdout.flush()


def _one(task):
    plt_cn, model_name = task
    sel = _G["sel"]
    trees = _G["trees"]
    conds = _G["conds"]
    row = sel.loc[plt_cn]

    t1 = trees.loc[trees["PLT_CN"] == str(plt_cn)].copy()
    if len(t1) < 3:
        return None

    cond_row = None
    c1 = conds.loc[conds["PLT_CN"].astype("int64") == plt_cn]
    if len(c1):
        cond_row = c1.iloc[0].to_dict()

    sid = f"NE_{plt_cn % 10_000_000:07d}"
    stand_df = R.build_standinit(sid, START_YEAR, VARIANT, cond_row,
                                 ecoregion=_G["ecomap"].get(plt_cn))
    tree_df = R.build_treeinit(t1, sid)
    if len(tree_df) == 0:
        return None
    tree_df = tree_df.drop(columns=["_cn"], errors="ignore")

    mdl = R.STANDALONE_MODELS[model_name]
    state_abbrev = R.STATECD_TO_ABBREV.get(int(row.get("STATECD", 0) or 0))
    restrict_dg = mdl.get("states_restrict_dg")
    if restrict_dg and state_abbrev not in restrict_dg:
        eff_extra_env = {k: v for k, v in mdl["extra_env"].items()
                         if k not in R._DG_KEYS}
        base_cc = mdl["calib_components"] or {}
        eff_calib = {**base_cc, "dg": True}
        dg_native = 0
    else:
        eff_extra_env = None
        eff_calib = None
        dg_native = 1

    df, _tl, err = R.run_stand(sid, stand_df, tree_df, VARIANT, model_name,
                               N_CYC, extra_env_override=eff_extra_env,
                               calib_override=eff_calib,
                               htg_kw="** BAKUZIS: no HTGMULT")
    if err or df is None or len(df) == 0:
        return None
    df = df.copy()
    df["PLOT"] = plt_cn
    df["arm"] = model_name
    df["site_cls"] = int(row["site_cls"])
    df["dens_cls"] = str(row["dens_cls"])
    df["SICOND"] = float(row["SICOND"])
    df["AGE0"] = float(row["STDAGE"])
    df["TPA0"] = float(row["TPA0"])
    df["dg_native"] = dg_native
    return df


def run_projection(sel_path: str, models: list[str], nproc: int,
                   outdir: str) -> None:
    states = R.VARIANTS[VARIANT]["states"]
    sel = pd.read_csv(sel_path)
    tasks = [(int(p), m) for m in models for p in sel["PLT_CN"]]
    print(f"[run] {len(sel)} plots x {len(models)} arms = {len(tasks)} runs, "
          f"{N_CYC} cycles x {CYCLE_LEN} yr = {N_CYC*CYCLE_LEN} yr horizon")
    sys.stdout.flush()

    _load_shared(sel_path, states)

    frames = []
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(_one, tasks, chunksize=4), 1):
            if r is not None:
                frames.append(r)
            if i % 200 == 0:
                print(f"  [{i}/{len(tasks)}] ok={len(frames)}")
                sys.stdout.flush()

    out = pd.concat(frames, ignore_index=True)
    p = os.path.join(outdir, "bakuzis_traj_raw.csv")
    out.to_csv(p, index=False)
    print(f"[run] wrote {len(out)} summary rows -> {p}")
    print(out.groupby("arm").size())
    print("columns:", list(out.columns))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["select", "run"])
    ap.add_argument("--plot-list", default="plot_list_3321.csv")
    ap.add_argument("--sel", default="bakuzis_plots.csv")
    ap.add_argument("--models", default="fvs_base,fvs_regional,organon,conus_spdep,conus_climate")
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--fia-dir", default=R.FIA_DIR)
    a = ap.parse_args()
    R.FIA_DIR = a.fia_dir
    os.makedirs(a.outdir, exist_ok=True)
    if a.stage == "select":
        select_plots(a.plot_list, a.sel)
    else:
        run_projection(a.sel, [m.strip() for m in a.models.split(",")],
                       a.nproc, a.outdir)


if __name__ == "__main__":
    main()
