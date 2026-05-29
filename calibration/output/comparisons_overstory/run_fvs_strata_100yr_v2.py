#!/usr/bin/env python3
"""run_fvs_strata_100yr_v2.py
v2: refined site index. Default SI=42 (balsam fir) was too optimistic
for these mature spruce/fir/mixed stands. Empirically observed dominant
heights of 25-35 ft on the SILC CFI at ~30-50 yr stand age put true SI
at ~35 ft. Lowering to SI=35 gives less optimistic growth projections.
"""
from __future__ import annotations
import os, sys, sqlite3, subprocess, tempfile
import pandas as pd, numpy as np

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", "/users/PUOM0008/crsfaaron/fvs-modern")
FVS_LIB_DIR  = os.environ.get("FVS_LIB_DIR", os.path.join(PROJECT_ROOT, "lib"))
sys.path.insert(0, "/users/PUOM0008/crsfaaron/silc_cfi")
from run_fvs_on_cfi import KEYFILE, calibrated_keywords

SITE_INDEX_REFINED = 35  # was 42; empirical estimate for mature spruce/fir

# Override build_standinit to use refined SI
def build_standinit_v2(sid, inv_year, variant):
    return pd.DataFrame([{
        "stand_id": sid, "variant": variant.upper(), "inv_year": int(inv_year),
        "latitude": 46.4628, "longitude": -68.4253, "region": 9,
        "forest": 0, "district": 0,
        "basal_area_factor": 0.0, "inv_plot_size": 1.0, "brk_dbh": 999.0,
        "num_plots": 1,
        "age": 60, "aspect": 0, "slope": 5, "elevft": 1000,
        "site_species": 12, "site_index": SITE_INDEX_REFINED,
        "state": 23, "county": 21,
        "forest_type": 121, "sam_wt": 1.0,
    }])

AGY_TO_SPCD = {
    "BF":  12, "RS":  97, "PB": 375, "YB": 371, "RM": 316,
    "SM": 318, "WA": 541, "RO": 833, "BS":  95, "WS":  91,
    "JP": 105, "RP": 125, "WP": 129, "EH": 261, "WC": 241,
    "BC": 531, "QA": 746, "GB": 934, "BA": 543, "AB": 531,
    "ST": 951, "TA": 547, "OH": 998, "OS": 299, "HW": 998,
    "PC": 117, "PR": 125, "AS": 543, "BT": 951, "EC": 241,
    "HH": 701, "NS": 94, "RB": 372, "RN": 94, "BP": 98,
    "SW": 998, "SB":  95, "TL": 547, "WB": 375,
}

def build_treeinit_strata(td, sid):
    rows = []
    for i, r in td.iterrows():
        spcd = AGY_TO_SPCD.get(str(r["Species"]), 998)
        dbh = float(r["DBH"])
        if not (np.isfinite(dbh) and dbh > 0): continue
        ht = float(r["Ht"]) if pd.notna(r["Ht"]) and r["Ht"] > 0 else 0.0
        tpa = float(r["TPA"])
        rows.append({"stand_id": sid, "plot_id": 1, "tree_id": i + 1,
                     "tree_count": round(tpa, 6), "species": spcd,
                     "diameter": round(dbh, 3), "ht": round(ht, 1), "crratio": 40})
    return pd.DataFrame(rows)

def run_one_100yr(sid, tree_df, inv_year, variant, config):
    binary = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}")
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, "FVS_Data.db")
        con = sqlite3.connect(db)
        build_standinit_v2(sid, inv_year, variant).to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        calib = calibrated_keywords(variant) if config == "calibrated" else "** DEFAULT"
        key = os.path.join(d, "strata.key")
        open(key, "w").write(KEYFILE.format(sid=sid, db=db, clen=5, ncyc=20, calib=calib))
        try:
            subprocess.run([binary, f"--keywordfile={key}"], cwd=d, capture_output=True, timeout=600)
        except subprocess.TimeoutExpired:
            return None
        try:
            df = pd.read_sql_query("SELECT * FROM FVS_Summary2", sqlite3.connect(db))
            df["variant"] = variant.upper(); df["config"]  = config
            return df
        except Exception:
            return None

def main():
    wd = os.getcwd()
    gr = pd.read_csv(os.path.join(wd, "GrownDB_byStrata_ALL.csv"))
    gr2023 = gr[(gr.Year == 2023) & gr.DBH.notna() & gr.TPA.notna() & (gr.TPA > 0)]
    stand_ids = sorted(gr2023["StandID"].unique())
    print(f"Running FVS-NE 100-yr v2 (SI={SITE_INDEX_REFINED}) on {len(stand_ids)} byStrata stands")
    rows = []
    for sid in stand_ids:
        td = gr2023[gr2023.StandID == sid]
        tree_df = build_treeinit_strata(td, sid)
        for variant in ("ne",):
            for config in ("default", "calibrated"):
                df = run_one_100yr(sid, tree_df, 2023, variant, config)
                if df is None or len(df) == 0:
                    print(f"  {sid} {variant} {config}: FAILED"); continue
                df["StandID"] = sid
                rows.append(df)
                last = df.iloc[-1]
                print(f"  {sid} {variant} {config}: yr2123 BA={last.get('BA',0):.1f} Cords={last.get('MCuFt',0)/79:.1f}")
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    out.to_csv(os.path.join(wd, "silc_strata_100yr_fvsne_si35_results.csv"), index=False)
    print(f"\nWrote {len(out)} rows for {out['StandID'].nunique() if len(out) else 0} stands")

if __name__ == "__main__":
    main()
