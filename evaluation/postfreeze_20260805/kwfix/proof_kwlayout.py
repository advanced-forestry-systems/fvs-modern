#!/usr/bin/env python3
"""Empirical proof of the BAMAX / BAIMULT / HTGMULT keyword column layouts.

Runs the FROZEN FVSne binary read-only on a synthetic stand with FVS keyword
echo enabled, emitting each keyword family in both the legacy 16-column form
written by config/config_loader.py and the corrected 10-column form, and reads
back what FVS itself echoed into the listing.  Nothing is asserted; every claim
comes from the listing.

Source facts established first, from src-converted/vbase/initre.f90:
  option 58 BAIMULT (activity 91): IDT=ARRAY(1) date, SPDECD(2,...) species,
      ARRAY(3) multiplier -> SCHEDULED activity, field 1 is the date.
  option 62 HTGMULT (activity 92): "GOTO 6005", identical processing to BAIMULT.
  option 66 BAMAX  : IF(ARRAY(1).GT.0.0) BAMAX=ARRAY(1) -> NOT scheduled,
      no date field and no species field at all.
"""
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

import pandas as pd

BINARY = "/users/PUOM0008/crsfaaron/fvs-modern/bin-balfix_20260804/FVSne"
OUTDIR = "/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/kwfix"


# ---------------------------------------------------------------- keyword forms
def old_bamax_scalar(v):
    return f"BAMAX           {v:10.1f}"


def old_bamax_species(i, v):
    return f"BAMAX           {i:10d}{v:10.1f}"


def new_bamax(v):
    return f"{'BAMAX':<10}{v:10.1f}"


def old_baimult(i, m):
    return f"BAIMULT         {i:10d}{m:10.4f}"


def new_baimult(i, m):
    return f"{'BAIMULT':<10}{'':10}{i:10d}{m:10.4f}"


def old_htgmult(i, m):
    return f"HTGMULT         {i:10d}{m:10.4f}"


def new_htgmult(i, m):
    return f"{'HTGMULT':<10}{'':10}{i:10d}{m:10.4f}"


def fields(rec):
    r = rec.ljust(70)
    return {
        "kw": r[0:8],
        "f1_date": r[10:20],
        "f2_species": r[20:30],
        "f3_value": r[30:40],
    }


# ---------------------------------------------------------------- synthetic stand
def standinit(sid):
    return pd.DataFrame([{
        "Stand_CN": sid, "StandPlot_CN": sid, "Stand_ID": sid, "StandPlot_ID": sid,
        "Variant": "NE", "Inv_Year": 2000, "Groups": "All_Stands", "AddFiles": "",
        "FVSKeywords": "", "Latitude": 46.5, "Longitude": -68.7, "Region": 9,
        "Forest": 12, "District": 1, "Compass": 0, "Basal_Area_Factor": -1.0,
        "Inv_Plot_Size": 1.0, "Brk_DBH": 5.0, "Num_Plots": 1, "NonStk_Plots": 0,
        "Sam_Wt": 1.0, "Stk_Pcnt": 100.0, "DG_Trans": 0, "DG_Measure": 0,
        "HTG_Trans": 0, "HTG_Measure": 0, "Mort_Measure": 0, "Max_BA": 0,
        "Max_SDI": 0, "Site_Species": "RS", "Site_Index": 42, "Model_Type": 4,
        "Physio_Region": 0, "Forest_Type": 0, "State": 23, "County": 1,
        "Fuel_Model": 0, "Slope": 5, "Aspect": 0, "Elevation": 10,
        "Age": 60, "PV_Code": "", "PV_Ref_Code": 0,
    }])


def treeinit(sid):
    rows = []
    for j in range(20):
        rows.append({
            "Stand_CN": sid, "StandPlot_CN": sid, "Stand_ID": sid,
            "StandPlot_ID": sid, "Plot_ID": 1, "Tree_ID": j + 1,
            "Tree_Count": 10.0, "History": 1, "Species": "RS",
            "DBH": 9.0 + 0.1 * j, "DG": "", "Ht": 50.0, "HtG": "", "HtTopK": "",
            "CrRatio": 5, "Damage1": 0, "Severity1": 0, "Damage2": 0,
            "Severity2": 0, "Damage3": 0, "Severity3": 0, "TreeValue": 1,
            "Prescription": 0, "Age": 60, "Slope": 5, "Aspect": 0,
            "PV_Code": "", "TopoCode": 0, "SitePrep": 0,
        })
    return pd.DataFrame(rows)


def run(sid, kw_block, keep_dir):
    tmp = tempfile.mkdtemp(prefix="kwlayout_")
    db = os.path.join(tmp, "FVS_Data.db")
    con = sqlite3.connect(db)
    standinit(sid).to_sql("fvs_standinit", con, if_exists="replace", index=False)
    treeinit(sid).to_sql("fvs_treeinit", con, if_exists="replace", index=False)
    con.commit()
    con.close()
    key = os.path.join(tmp, "run.key")
    with open(key, "w") as fh:
        fh.write(
            "STDIDENT\n" + sid + "\n"
            "DATABASE\nDSNIN\n" + db + "\nDSNOUT\n" + db + "\n"
            "STANDSQL\nSELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'\nENDSQL\n"
            "TREESQL\nSELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'\nENDSQL\n"
            "END\n"
            "DATABASE\nSUMMARY            2\nEND\n"
            "FIAVBC\n"
            "TIMEINT            0        10\n"
            "NUMCYCLE           2\n"
            + kw_block + "\n"
            "PROCESS\nSTOP\n"
        )
    p = subprocess.run([BINARY, f"--keywordfile={key}"], cwd=tmp,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out = ""
    for cand in os.listdir(tmp):
        if cand.endswith(".out"):
            out = open(os.path.join(tmp, cand), errors="replace").read()
            break
    os.makedirs(keep_dir, exist_ok=True)
    shutil.copy(key, os.path.join(keep_dir, "run.key"))
    with open(os.path.join(keep_dir, "run.out"), "w") as fh:
        fh.write(out)
    shutil.rmtree(tmp, ignore_errors=True)
    return p.returncode, out


def echo_lines(out, tokens):
    keep = []
    for i, ln in enumerate(out.splitlines()):
        if any(t in ln for t in tokens):
            keep.append((i + 1, ln.rstrip()))
    return keep


def main():
    res = {}
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("FVS keyword column-layout proof, frozen binary")
    say("binary : " + BINARY)
    say("")

    # --- BAMAX: legacy scalar, legacy per-species, then corrected -------------
    kb = "\n".join([
        old_bamax_scalar(250.0),
        old_bamax_species(2, 250.0),
        new_bamax(250.0),
    ])
    say("=== BAMAX keyword block written ===")
    for r in kb.splitlines():
        say("  |" + r + "|   " + json.dumps(fields(r)))
    rc, out = run("KWBAMAX", kb, os.path.join(OUTDIR, "case_bamax"))
    say(f"returncode = {rc}")
    say("--- listing echo of MAXIMUM BASAL AREA ---")
    ech = echo_lines(out, ["MAXIMUM BASAL AREA"])
    for n, ln in ech:
        say(f"  line {n}: {ln}")
    vals = []
    for _, ln in ech:
        try:
            vals.append(float(ln.split("=")[-1]))
        except ValueError:
            vals.append(None)
    res["bamax_echoed_values_in_order"] = vals
    res["bamax_legacy_scalar_bound"] = bool(vals and vals[0] not in (None, 0.0))
    res["bamax_legacy_species_bound"] = bool(len(vals) > 1 and vals[1] not in (None, 0.0))
    res["bamax_corrected_bound"] = bool(len(vals) > 2 and abs((vals[2] or 0) - 250.0) < 0.01)
    say("")

    # --- BAIMULT / HTGMULT ---------------------------------------------------
    kb = "\n".join([
        old_baimult(2, 1.8914),
        new_baimult(2, 1.8914),
        old_htgmult(2, 1.8914),
        new_htgmult(2, 1.8914),
        old_baimult(12, 12.3456),
        new_baimult(12, 12.3456),
    ])
    say("=== BAIMULT / HTGMULT keyword block written ===")
    for r in kb.splitlines():
        say("  |" + r + "|   " + json.dumps(fields(r)))
    rc, out = run("KWMULT", kb, os.path.join(OUTDIR, "case_mult"))
    say(f"returncode = {rc}")
    say("--- listing echo of BAIMULT / HTGMULT / FVS errors ---")
    ech = echo_lines(out, ["BAIMULT", "HTGMULT", "FVS0", "FVS1"])
    for n, ln in ech:
        say(f"  line {n}: {ln}")
    res["mult_echo_lines"] = [ln for _, ln in ech]
    res["mult_returncode"] = rc
    res["mult_fvs_errors"] = [ln for _, ln in ech if "ERROR" in ln]
    say("")

    with open(os.path.join(OUTDIR, "PROOF_kwlayout.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUTDIR, "proof_kwlayout.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote PROOF_kwlayout.txt and proof_kwlayout.json")


if __name__ == "__main__":
    main()
