#!/usr/bin/env python3
"""
run_fvs_treelist_v2.py
======================
Per-tree FVS output by parsing the TREELIST .trl text file.
"""
from __future__ import annotations
import os, re, sys, sqlite3, subprocess, tempfile
import pandas as pd

sys.path.insert(0, "/users/PUOM0008/crsfaaron/silc_cfi")
from run_fvs_on_cfi import build_standinit, build_treeinit_cfi, FVS_LIB_DIR

def build_keyfile(sid, db, ncyc, clen, calib, inv_year):
    # Schedule TREELIST at each cycle year so the final-cycle tree list is emitted
    tl_lines = "\n".join(f"TREELIST           {inv_year + k}" for k in range(1, ncyc + 1))
    return f"""STDIDENT
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
{tl_lines}
DATABASE
SUMMARY            2
END
TIMEINT            0         {clen}
NUMCYCLE          {ncyc}
{calib}
PROCESS
STOP
"""

# Section header pattern: "COMPLETE TREE LIST -- STAND: <sid>" with "YEAR: <yyyy>"
HDR_RE = re.compile(r"COMPLETE TREE LIST.*STAND:\s*(\S+).*YEAR:\s*(\d+)")

def parse_trl(path):
    """Yields (year, df) for each cycle section in the .trl file."""
    if not os.path.exists(path):
        return
    with open(path) as f:
        lines = f.readlines()
    # locate section headers and the data block following each
    sections = []
    for i, ln in enumerate(lines):
        m = HDR_RE.search(ln)
        if m:
            sections.append((int(m.group(2)), i))
    sections.append((None, len(lines)))
    for k in range(len(sections) - 1):
        year, start = sections[k]
        end = sections[k+1][1]
        block = lines[start:end]
        rows = []
        for ln in block:
            # data row pattern: starts with tree number then values
            # Use fixed-column layout sniffed from header positions
            if not ln.strip() or ln.startswith("-") or ln.startswith("="):
                continue
            # Skip header rows (contain TREE NUMBER or letters)
            if "TREE" in ln and "NUMBER" in ln:
                continue
            if "INDX" in ln or "SP CD" in ln:
                continue
            # Data row: try to parse the species code (2-3 chars after tree numbers)
            parts = ln.split()
            if len(parts) < 8:
                continue
            try:
                tree_num = int(parts[0])
            except ValueError:
                continue
            # Position of SP CD column varies; use it as third token typically
            sp_code = parts[2]
            if not sp_code.isalpha():
                continue
            try:
                tpa  = float(parts[6])
                dbh  = float(parts[8])
                ht   = float(parts[10])
            except (ValueError, IndexError):
                continue
            rows.append((sp_code, tpa, dbh, ht))
        if rows:
            df = pd.DataFrame(rows, columns=["SP","TPA","DBH_in","HT_ft"])
            df["Year"] = year
            yield year, df

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
        build_standinit(sid, inv_year, variant).to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        calib = (calibrated_keywords(variant) if config == "calibrated" else "** DEFAULT")
        key = os.path.join(d, "cfi.key")
        open(key, "w").write(build_keyfile(sid, db, ncyc, clen, calib, inv_year))
        try:
            subprocess.run([binary, f"--keywordfile={key}"], cwd=d, capture_output=True, timeout=180)
        except subprocess.TimeoutExpired:
            return None
        trl = os.path.join(d, "cfi.trl")
        if not os.path.exists(trl):
            return None
        # Return only the last cycle's per-tree records (year = inv_year + ncyc)
        target_year = inv_year + ncyc
        last_df = None
        for yr, df in parse_trl(trl):
            if yr == target_year:
                last_df = df
        return last_df

def main():
    wd = os.getcwd()
    manifest = pd.read_csv(os.path.join(wd, "silc_cfi_pair_summary.csv"))
    rows = []
    print(f"[fvs tl v2] running {len(manifest)} pairs x 2 variants x 2 configs")
    for i, pr in manifest.iterrows():
        tf = os.path.join(wd, pr["tree_list_file"])
        if not os.path.exists(tf): continue
        td = pd.read_csv(tf)
        if len(td) == 0: continue
        sid = f"CFI_{int(pr['PLOT']):04d}_{int(pr['YEAR_PREV'])}"
        tree_df = build_treeinit_cfi(td, sid)
        if len(tree_df) == 0: continue
        ncyc = int(pr["PERIOD_YR"])
        for variant in ("ne", "acd"):
            for config in ("default", "calibrated"):
                df = run_one(sid, tree_df, int(pr["YEAR_PREV"]),
                             variant, config, ncyc, clen=1)
                if df is None or len(df) == 0:
                    continue
                for _, r in df.iterrows():
                    rows.append({
                        "PLOT": int(pr["PLOT"]),
                        "YEAR_PREV": int(pr["YEAR_PREV"]),
                        "YEAR_CURR": int(pr["YEAR_CURR"]),
                        "variant": variant.upper(),
                        "config":  config,
                        "SP":      r["SP"],
                        "DBH_in":  float(r["DBH_in"]),
                        "HT_ft":   float(r["HT_ft"]),
                        "EXPF_ac": float(r["TPA"]),
                    })
        if i % 4 == 0:
            print(f"  done pair {i+1}/{len(manifest)}")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(wd, "silc_cfi_fvs_treelist.csv"), index=False)
    print(f"[fvs tl v2] wrote silc_cfi_fvs_treelist.csv: {len(out)} per-tree rows")

if __name__ == "__main__":
    main()
