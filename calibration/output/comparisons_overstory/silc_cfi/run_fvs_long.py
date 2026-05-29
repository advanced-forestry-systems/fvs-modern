#!/usr/bin/env python3
"""run_fvs_long.py - Cardinal driver for long-horizon FVS-NE/ACD"""
from __future__ import annotations
import os, sys, sqlite3, subprocess, tempfile
import pandas as pd, numpy as np
sys.path.insert(0, "/users/PUOM0008/crsfaaron/silc_cfi")
from run_fvs_on_cfi import (build_standinit, build_treeinit_cfi, KEYFILE,
                             calibrated_keywords, FVS_LIB_DIR)

def run_one(sid, tree_df, inv_year, variant, config, ncyc, clen=1):
    binary = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}")
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, "FVS_Data.db")
        con = sqlite3.connect(db)
        build_standinit(sid, inv_year, variant).to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        calib = calibrated_keywords(variant) if config == "calibrated" else "** DEFAULT"
        key = os.path.join(d, "cfi.key")
        open(key, "w").write(KEYFILE.format(sid=sid, db=db, clen=clen, ncyc=ncyc, calib=calib))
        try:
            subprocess.run([binary, f"--keywordfile={key}"], cwd=d, capture_output=True, timeout=300)
        except subprocess.TimeoutExpired:
            return None
        try:
            df = pd.read_sql_query("SELECT * FROM FVS_Summary2", sqlite3.connect(db))
            df["variant"] = variant.upper(); df["config"] = config
            return df
        except Exception:
            return None

def main():
    wd = os.getcwd()
    manifest = pd.read_csv(os.path.join(wd, "silc_cfi_longhorizon_pairs.csv"))
    print(f"[fvs long] {len(manifest)} pairs x 2 variants x 2 configs")
    rows = []
    for i, pr in manifest.iterrows():
        tf = os.path.join(wd, pr["tree_list_file"])
        if not os.path.exists(tf): continue
        td = pd.read_csv(tf)
        if len(td) == 0: continue
        sid = f"CFI_{int(pr['PLOT']):04d}_{int(pr['YEAR_PREV'])}"
        tree_df = build_treeinit_cfi(td, sid)
        ncyc = int(pr["PERIOD_YR"])
        for variant in ("ne", "acd"):
            for config in ("default", "calibrated"):
                df = run_one(sid, tree_df, int(pr["YEAR_PREV"]),
                             variant, config, ncyc, clen=1)
                if df is None or len(df) == 0: continue
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
        print(f"  pair {i+1}/{len(manifest)} plot {int(pr['PLOT'])} ({int(pr['PERIOD_YR'])}yr): done")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(wd, "silc_cfi_long_fvs_results.csv"), index=False)
    print(f"[fvs long] wrote {len(out)} rows")

if __name__ == "__main__":
    main()
