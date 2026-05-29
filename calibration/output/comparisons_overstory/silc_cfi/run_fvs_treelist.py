#!/usr/bin/env python3
"""
run_fvs_treelist.py
===================
Patched FVS driver that adds TREELIST output to capture per-tree
predicted year_curr DBH, HT, EXPF, and species. Writes one CSV per
configuration. The standalone FVS binary writes FVS_TreeList table
(or FVS_TreeList_East depending on variant) via TREELIDB keyword;
we pull the final cycle's per-tree records.

Output: silc_cfi_fvs_treelist.csv (long form: PLOT, YEAR_PREV,
YEAR_CURR, variant, config, SPCD, DBH_in, HT_ft, EXPF_ac).
"""
from __future__ import annotations
import os, sys, sqlite3, subprocess, tempfile
import pandas as pd
import numpy as np

sys.path.insert(0, "/users/PUOM0008/crsfaaron/silc_cfi")
from run_fvs_on_cfi import build_standinit, build_treeinit_cfi, FVS_LIB_DIR

# Keyfile with TREELIDB enabled to dump per-tree records
KEYFILE_TL = """STDIDENT
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
TREELIDB           2         2
END
TIMEINT            0         {clen}
NUMCYCLE          {ncyc}
{calib}
PROCESS
STOP
"""

def calibrated_keywords(variant):
    try:
        sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-modern")
        from config.config_loader import FvsConfigLoader
        return FvsConfigLoader(variant.lower(), version="calibrated",
                               config_dir="/users/PUOM0008/crsfaaron/fvs-modern/config"
                               ).generate_keywords(include_comments=False)
    except Exception:
        return "** DEFAULT"

def run_one(sid, tree_df, inv_year, variant, config, ncyc, clen=1):
    binary = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}")
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, "FVS_Data.db")
        con = sqlite3.connect(db)
        build_standinit(sid, inv_year, variant).to_sql(
            "fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        calib = (calibrated_keywords(variant) if config == "calibrated"
                 else "** DEFAULT")
        key = os.path.join(d, "cfi.key")
        open(key, "w").write(KEYFILE_TL.format(
            sid=sid, db=db, clen=clen, ncyc=ncyc, calib=calib))
        try:
            subprocess.run([binary, f"--keywordfile={key}"], cwd=d,
                           capture_output=True, timeout=120)
        except subprocess.TimeoutExpired:
            return None
        con = sqlite3.connect(db)
        tables = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
        # Find the FVS_TreeList table (may have variant suffix)
        tl_table = next((t for t in tables if "TreeList" in t), None)
        if not tl_table:
            return None
        try:
            df = pd.read_sql(f"SELECT * FROM {tl_table}", con)
            con.close()
            return df
        except Exception:
            return None

def main():
    wd = os.getcwd()
    manifest = pd.read_csv(os.path.join(wd, "silc_cfi_pair_summary.csv"))
    rows = []
    print(f"[fvs tl] running {len(manifest)} pairs x 2 variants x 2 configs")
    for i, pr in manifest.iterrows():
        tf = os.path.join(wd, pr["tree_list_file"])
        if not os.path.exists(tf):
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
                df = run_one(sid, tree_df, int(pr["YEAR_PREV"]),
                             variant, config, ncyc, clen=1)
                if df is None or len(df) == 0:
                    continue
                # Keep only final-cycle trees (max Year)
                last_year = df["Year"].max() if "Year" in df.columns else None
                if last_year is not None:
                    df = df[df["Year"] == last_year]
                # Common FVS_TreeList columns: SpeciesPLANTS, SpeciesFIA, DBH, Ht, TPA
                # Pick whichever column names exist
                def col(*names):
                    for n in names:
                        if n in df.columns:
                            return n
                    return None
                c_sp  = col("SpeciesFIA","FIA_CODE","SPCD","Species")
                c_d   = col("DBH","Dbh","DIA")
                c_h   = col("Ht","HT","Height")
                c_e   = col("TPA","Tpa","TreeCount","Stems")
                for _, r in df.iterrows():
                    rows.append({
                        "PLOT": int(pr["PLOT"]),
                        "YEAR_PREV": int(pr["YEAR_PREV"]),
                        "YEAR_CURR": int(pr["YEAR_CURR"]),
                        "variant":  variant.upper(),
                        "config":   config,
                        "SPCD":    int(r[c_sp]) if c_sp and pd.notna(r[c_sp]) else None,
                        "DBH_in": float(r[c_d]) if c_d and pd.notna(r[c_d]) else None,
                        "HT_ft":  float(r[c_h]) if c_h and pd.notna(r[c_h]) else None,
                        "EXPF_ac": float(r[c_e]) if c_e and pd.notna(r[c_e]) else None,
                    })
        if i % 4 == 0:
            print(f"  done pair {i+1}/{len(manifest)}")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(wd, "silc_cfi_fvs_treelist.csv"), index=False)
    print(f"[fvs tl] wrote silc_cfi_fvs_treelist.csv: {len(out)} per-tree rows")

if __name__ == "__main__":
    main()
