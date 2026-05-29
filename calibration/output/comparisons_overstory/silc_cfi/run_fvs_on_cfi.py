#!/usr/bin/env python3
"""
run_fvs_on_cfi.py
=================
Cardinal-side driver. Runs FVS-NE (Northeast) and FVS-ACD (Acadian)
in BOTH default and calibrated configurations on the 24 SILC CFI
pair tree lists.

Pattern follows the working magplot_fvs_runner.py: standalone FVS
binary, DATABASE keyfile, one stand per subprocess. Adapted for CFI:
imperial input units, FIA SPCDs already present, PERIOD_YR cycles at
clen=1 yr per cycle.

Inputs (in current working directory):
  silc_cfi_pair_summary.csv         24 pair manifest
  pair_input/pair_<P>_<Yprev>_tree.csv  per-pair tree lists

Output:
  silc_cfi_fvs_results.csv  one row per (pair, variant, config) with
                            BA, TPA, QMD, MCuFt, BdFt at year_curr
"""
from __future__ import annotations
import os, sys, sqlite3, subprocess, tempfile
import pandas as pd
import numpy as np

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", "/users/PUOM0008/crsfaaron/fvs-modern")
FVS_LIB_DIR  = os.environ.get("FVS_LIB_DIR", os.path.join(PROJECT_ROOT, "lib"))
CONFIG_DIR   = os.environ.get("FVS_CONFIG_DIR", os.path.join(PROJECT_ROOT, "config"))
ACRES_PER_HA = 2.4710538147

KEYFILE = """STDIDENT
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
TIMEINT            0         {clen}
NUMCYCLE          {ncyc}
{calib}
PROCESS
STOP
"""

def build_standinit(sid, inv_year, variant):
    return pd.DataFrame([{
        "stand_id": sid, "variant": variant.upper(), "inv_year": int(inv_year),
        "latitude": 46.5, "longitude": -68.7, "region": 9,
        "forest": 0, "district": 0,
        "basal_area_factor": 0.0, "inv_plot_size": 0.2, "brk_dbh": 4.5,
        "num_plots": 1,
        "age": 60, "aspect": 0, "slope": 5, "elevft": 1000,
        "site_species": 12, "site_index": 42, "state": 23, "county": 21,
        "forest_type": 121, "sam_wt": 1.0,
    }])

def build_treeinit_cfi(td, sid):
    """CFI pair_input tree list is already imperial: DIA_IN, HT_FT, EXPF=5/ac.
    Build FVS treeinit rows."""
    rows = []
    for i, r in td.iterrows():
        spcd = int(r["SPCD"]) if pd.notna(r["SPCD"]) else 12
        dbh  = float(r["DBH"])
        if not (np.isfinite(dbh) and dbh > 0):
            continue
        ht   = float(r["HT"]) if pd.notna(r["HT"]) and r["HT"] > 0 else 0.0
        # CFI EXPF is 5 trees/ac per tree (1/5 ac plot)
        tpa  = float(r["EXPF"])
        rows.append({
            "stand_id": sid, "plot_id": 1, "tree_id": int(r["TREE"]),
            "tree_count": round(tpa, 5), "species": spcd,
            "diameter": round(dbh, 3), "ht": round(ht, 1),
            "crratio": 40,
        })
    return pd.DataFrame(rows)

def calibrated_keywords(variant):
    try:
        sys.path.insert(0, PROJECT_ROOT)
        from config.config_loader import FvsConfigLoader
        return FvsConfigLoader(variant.lower(), version="calibrated",
                               config_dir=CONFIG_DIR).generate_keywords(include_comments=False)
    except Exception as e:
        sys.stderr.write(f"  calibrated kw unavailable for {variant}: {e}\n")
        return "** DEFAULT (calibrated config not found)"

def run_stand(sid, tree_df, inv_year, variant, config, ncyc, clen=1):
    binary = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}")
    if not os.path.exists(binary):
        return None, f"binary not found: {binary}"
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, "FVS_Data.db")
        con = sqlite3.connect(db)
        build_standinit(sid, inv_year, variant).to_sql("fvs_standinit",
                                                        con, if_exists="replace",
                                                        index=False)
        tree_df.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        calib = (calibrated_keywords(variant) if config == "calibrated"
                 else "** DEFAULT PARAMETERS")
        key = os.path.join(d, "cfi.key")
        with open(key, "w") as f:
            f.write(KEYFILE.format(sid=sid, db=db, clen=clen, ncyc=ncyc, calib=calib))
        try:
            subprocess.run([binary, f"--keywordfile={key}"], cwd=d,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           timeout=120)
        except subprocess.TimeoutExpired:
            return None, "timeout"
        try:
            con = sqlite3.connect(db)
            df = pd.read_sql_query("SELECT * FROM FVS_Summary2", con)
            con.close()
            df["variant"] = variant.upper()
            df["config"] = config
            return df, None
        except Exception as e:
            return None, f"no FVS_Summary2: {e}"

def main():
    wd = os.getcwd()
    manifest = pd.read_csv(os.path.join(wd, "silc_cfi_pair_summary.csv"))
    print(f"[cfi] running {len(manifest)} pairs x 2 variants x 2 configs "
          f"= {len(manifest)*4} FVS runs")
    rows = []
    for i, pr in manifest.iterrows():
        tf = os.path.join(wd, pr["tree_list_file"])
        if not os.path.exists(tf):
            print(f"  pair {i+1}/{len(manifest)}: missing {pr['tree_list_file']}")
            continue
        td = pd.read_csv(tf)
        if len(td) == 0:
            continue
        sid = f"CFI_{int(pr['PLOT']):04d}_{int(pr['YEAR_PREV'])}"
        tree_df = build_treeinit_cfi(td, sid)
        if len(tree_df) == 0:
            continue
        ncyc = int(pr["PERIOD_YR"])
        for variant in ("ne", "acd"):
            for config in ("default", "calibrated"):
                df, err = run_stand(sid, tree_df, int(pr["YEAR_PREV"]),
                                    variant, config, ncyc, clen=1)
                if err:
                    print(f"  {sid} {variant}/{config}: ERROR {err}")
                    continue
                if len(df) == 0:
                    print(f"  {sid} {variant}/{config}: empty Summary2")
                    continue
                # Last cycle = year_curr
                yn = df.iloc[-1]
                rows.append({
                    "PLOT": int(pr["PLOT"]),
                    "YEAR_PREV": int(pr["YEAR_PREV"]),
                    "YEAR_CURR": int(pr["YEAR_CURR"]),
                    "PERIOD_YR": ncyc,
                    "variant": variant.upper(),
                    "config":  config,
                    "BA_PRED_ft2ac": float(yn.get("BA", np.nan)),
                    "TPA_PRED":      float(yn.get("Tpa", np.nan)),
                    "QMD_PRED_in":   float(yn.get("QMD", np.nan)),
                    "MCuFt_PRED":    float(yn.get("MCuFt", np.nan)),
                    "BdFt_PRED":     float(yn.get("BdFt", np.nan)),
                    "BA_OBS_PREV":   float(pr["BA_PREV_FT2AC"]),
                    "BA_OBS_CURR":   float(pr["BA_CURR_FT2AC"]),
                })
        print(f"  pair {i+1}/{len(manifest)} {sid}: done")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(wd, "silc_cfi_fvs_results.csv"), index=False)
    print(f"\n[cfi] wrote silc_cfi_fvs_results.csv: {len(out)} rows")
    for (v, c), grp in out.groupby(["variant", "config"]):
        ok = grp["BA_PRED_ft2ac"].notna() & grp["BA_OBS_CURR"].notna()
        if ok.sum() > 0:
            bias_pct = 100*(grp.loc[ok, "BA_PRED_ft2ac"].mean()
                            / grp.loc[ok, "BA_OBS_CURR"].mean() - 1)
            rmse = float(np.sqrt(((grp.loc[ok, "BA_PRED_ft2ac"]
                                  - grp.loc[ok, "BA_OBS_CURR"])**2).mean()))
            print(f"  FVS-{v} {c}: n={ok.sum()}  BA bias {bias_pct:+.2f}%  "
                  f"RMSE {rmse:.2f} ft^2/ac")

if __name__ == "__main__":
    main()
