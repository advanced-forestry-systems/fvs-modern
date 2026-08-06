#!/usr/bin/env python3
"""
proof_maxcyc.py -- prove the MAXCYC / NUMCYCLE silent-truncation guard with a
real FVS run on a SYNTHETIC stand (no FIA plot data).

Frozen binary: /users/PUOM0008/crsfaaron/fvs-modern/bin-balfix_20260804/FVSne
(read-only, never overwritten).

Emits PROOF.txt and proof.json in WORK.
"""
from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

import pandas as pd

WORK = "/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/maxcyc"
DRIVER_DIR = "/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804"
TOOLS_DIR = "/users/PUOM0008/crsfaaron/fvs-modern/tools"
PRGPRM = "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/base/PRGPRM.f90"

sys.path.insert(0, TOOLS_DIR)
sys.path.insert(0, DRIVER_DIR)
import run_driver_final as RD           # noqa: E402  read-only import
import fvs_run_guards as G              # noqa: E402

BINARY = RD.VARIANT_BINS["ne"]
ARM = "fvs_base"

LOG = []


def say(s=""):
    print(s)
    LOG.append(str(s))


# ---------------------------------------------------------------------------
# Synthetic stand (construction copied from stress_harness.py)
# ---------------------------------------------------------------------------
PIRU = 97


def standinit(sid, site_index=42, elev=1000, age=60, slope=5, inv_year=2000):
    return pd.DataFrame([{
        "stand_id": sid, "variant": "NE", "inv_year": int(inv_year),
        "latitude": 46.5, "longitude": -68.7, "region": 9, "forest": 0,
        "district": 0, "basal_area_factor": 0.0, "inv_plot_size": 1.0,
        "brk_dbh": 999.0, "num_plots": 1, "age": int(age), "aspect": 0,
        "slope": int(slope), "elevft": int(elev), "site_species": 12,
        "site_index": float(site_index), "state": 23, "county": 0,
        "forest_type": 121, "sam_wt": 1.0, "ecoregion": "M210",
    }])


def treeinit(sid, recs):
    rows = []
    for i, (sp, dbh, tpa, ht, cr) in enumerate(recs, start=1):
        rows.append({"stand_id": sid, "plot_id": 1, "tree_id": i,
                     "tree_count": round(float(tpa), 4), "species": int(sp),
                     "diameter": round(float(dbh), 2), "ht": round(float(ht), 1),
                     "crratio": int(max(15, min(99, cr)))})
    return pd.DataFrame(rows)


def cycle_kw(clen, ncyc):
    """Annual convention: TIMEINT 0 <clen> / NUMCYCLE <ncyc>."""
    return f"TIMEINT            0{clen:>10d}\nNUMCYCLE  {ncyc:>10d}"


# ---------------------------------------------------------------------------
# Raw (unguarded) FVS invocation -- this is exactly what run_driver_final does,
# except we keep stdout/stderr and the listing instead of discarding them.
# ---------------------------------------------------------------------------
def raw_run(sid, ncyc, clen=1, keep_dir=None):
    mdl = RD.STANDALONE_MODELS[ARM]
    env = os.environ.copy()
    env.update(mdl["extra_env"])
    calib = RD.get_calib_kw("ne", mdl["calib_components"])

    tmp = tempfile.mkdtemp(prefix="maxcycproof_", dir="/tmp")
    db = os.path.join(tmp, "FVS_Data.db")
    con = sqlite3.connect(db)
    standinit(sid).to_sql("fvs_standinit", con, if_exists="replace", index=False)
    treeinit(sid, [(PIRU, 9.0, 200.0, 0.0, 50)]).to_sql(
        "fvs_treeinit", con, if_exists="replace", index=False)
    con.close()

    key = os.path.join(tmp, "run.key")
    with open(key, "w") as fh:
        fh.write(RD.KEYFILE.format(sid=sid, db=db,
                                   cycle_kw=cycle_kw(clen, ncyc),
                                   calib_kw=calib, htg_kw=""))

    p = subprocess.run([BINARY, f"--keywordfile={key}"], cwd=tmp,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       env=env, timeout=900)
    out = os.path.join(tmp, "run.out")
    rows = G.summary_row_count(db)
    res = {"rc": p.returncode, "key": key, "out": out, "db": db, "rows": rows,
           "tmp": tmp,
           "stderr": p.stderr.decode("utf8", "replace")[:400]}
    if keep_dir:
        os.makedirs(keep_dir, exist_ok=True)
        for f in ("run.key", "run.out"):
            src = os.path.join(tmp, f)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(keep_dir, f))
        res["key"] = os.path.join(keep_dir, "run.key")
        res["out"] = os.path.join(keep_dir, "run.out")
    return res


# ---------------------------------------------------------------------------
def main():
    proof = {}

    say("=" * 78)
    say("PROOF: MAXCYC / NUMCYCLE silent-truncation guard")
    say("date: 2026-08-05/06   binary: " + BINARY)
    say("stand: SYNTHETIC even-aged Picea rubens monoculture, 200 TPA, 9.0 in "
        "DBH, SI 42, no FIA data")
    say("=" * 78)

    # ---- 5. MAXCYC parse ---------------------------------------------------
    say("")
    say("[5] PARSED MAXCYC FROM SOURCE")
    say("-" * 78)
    maxcyc = G.read_maxcyc(PRGPRM)
    say("$ grep -n 'PARAMETER (MAXCYC' %s" % PRGPRM)
    with open(PRGPRM) as fh:
        for i, ln in enumerate(fh, 1):
            if "MAXCYC=" in ln.upper().replace(" ", "").replace("MAXCYC =", "MAXCYC="):
                say("%d:%s" % (i, ln.rstrip()))
    say("G.read_maxcyc(...) -> %d   (MAXCYC_DEFAULT fallback = %d)"
        % (maxcyc, G.MAXCYC_DEFAULT))
    say("MATCHES 40: %s" % (maxcyc == 40))
    proof["maxcyc_parsed"] = maxcyc

    # ---- 2. CONTROL: unguarded NUMCYCLE 41 ---------------------------------
    say("")
    say("[2] CONTROL -- UNGUARDED RUN, NUMCYCLE 41 (the silent failure)")
    say("-" * 78)
    keep41 = os.path.join(WORK, "case_numcycle41")
    r41 = raw_run("SYN41", 41, clen=1, keep_dir=keep41)
    say("$ %s --keywordfile=run.key    # TIMEINT 0 1 / NUMCYCLE 41" % BINARY)
    say("return code            : %s   (normal FVS STOP is 10)" % r41["rc"])
    say("FVS_Summary2 row count : %s" % r41["rows"])
    codes41 = G.scan_listing(r41["out"])
    say("FVS error codes in listing:")
    seen = []
    for c in codes41:
        if c["code"] not in [x["code"] for x in seen]:
            seen.append(c)
    for c in seen:
        say("   %s %-7s line %d: %s"
            % (c["code"], c["severity"], c["lineno"], c["line"].strip()))
    if not codes41:
        say("   (none found)")
    fvs04 = next((c for c in codes41 if c["code"] == "FVS04"), None)
    say("stderr head: %s" % (r41["stderr"].replace("\n", " ")[:200] or "(empty)"))
    say("")
    say("=> An unguarded driver that DEVNULLs stdout/stderr and ignores the "
        "return code sees a well-formed %s-row summary and reports success."
        % r41["rows"])
    proof["unguarded_returncode"] = r41["rc"]
    proof["unguarded_summary_rows"] = r41["rows"]
    proof["unguarded_fvs04_line"] = fvs04["line"].strip() if fvs04 else None
    proof["unguarded_all_codes"] = sorted({c["code"] for c in codes41})

    # ---- 3a. GUARDED PRE-SUBMIT --------------------------------------------
    say("")
    say("[3a] GUARDED PATH -- PRE-SUBMIT, FVS IS NEVER LAUNCHED")
    say("-" * 78)
    say(">>> fvs_run_guards.validate_numcycle(41)")
    presubmit = False
    try:
        G.validate_numcycle(41, maxcyc=maxcyc)
        say("NO ERROR RAISED -- GUARD FAILED")
    except G.FVSKeywordError as exc:
        presubmit = True
        say("FVSKeywordError: %s" % exc)
    say(">>> fvs_run_guards.validate_keyword_file(%s)" % r41["key"])
    presubmit_kf = False
    try:
        G.validate_keyword_file(r41["key"], maxcyc=maxcyc)
        say("NO ERROR RAISED -- GUARD FAILED")
    except G.FVSKeywordError as exc:
        presubmit_kf = True
        say("FVSKeywordError: %s" % str(exc)[:200] + " ...")
    proof["guarded_presubmit_rejected"] = bool(presubmit and presubmit_kf)

    # ---- 3b. GUARDED POST-RUN ----------------------------------------------
    say("")
    say("[3b] GUARDED PATH -- POST-RUN, ON THE CONTROL RUN'S ARTIFACTS")
    say("-" * 78)
    say(">>> fvs_run_guards.check_run(returncode=%s, out_path=..., "
        "expected_cycles=41, summary_rows=%s)" % (r41["rc"], r41["rows"]))
    postrun = False
    try:
        G.check_run(r41["rc"], r41["out"], expected_cycles=41,
                    summary_rows=r41["rows"], stderr_text=r41["stderr"])
        say("NO ERROR RAISED -- GUARD FAILED")
    except G.FVSRunError as exc:
        postrun = True
        say("FVSRunError: %s" % exc)
    proof["guarded_postrun_detected"] = bool(postrun)

    # ---- 4. BOUNDARY: NUMCYCLE 40 ------------------------------------------
    say("")
    say("[4] BOUNDARY -- NUMCYCLE 40 THROUGH THE GUARDED PATH (must PASS)")
    say("-" * 78)
    say(">>> fvs_run_guards.validate_numcycle(40)")
    try:
        G.validate_numcycle(40, maxcyc=maxcyc)
        say("pre-submit OK, launching FVS")
        pre40 = True
    except G.FVSKeywordError as exc:
        pre40 = False
        say("UNEXPECTED pre-submit reject: %s" % exc)

    keep40 = os.path.join(WORK, "case_numcycle40")
    r40 = raw_run("SYN40", 40, clen=1, keep_dir=keep40)
    say("$ %s --keywordfile=run.key    # TIMEINT 0 1 / NUMCYCLE 40" % BINARY)
    say("return code            : %s" % r40["rc"])
    say("FVS_Summary2 row count : %s   (expected 41 = 40 cycles + 1)" % r40["rows"])
    codes40 = G.scan_listing(r40["out"])
    say("FVS codes in listing: %s"
        % (sorted({(c["code"], c["severity"]) for c in codes40}) or "(none)"))
    for c in codes40[:5]:
        say("   %s %-7s line %d: %s"
            % (c["code"], c["severity"], c["lineno"], c["line"].strip()))
    say("NOTE: FVS03 is a WARNING (synthetic stand has no real forest code), "
        "not an ERROR.  check_run fails on ERROR severity only, so this "
        "legal run passes while still reporting the warning.")
    say(">>> fvs_run_guards.check_run(returncode=%s, out_path=..., "
        "expected_cycles=40, summary_rows=%s)" % (r40["rc"], r40["rows"]))
    b40 = False
    try:
        res = G.check_run(r40["rc"], r40["out"], expected_cycles=40,
                          summary_rows=r40["rows"])
        b40 = True
        say("POST-RUN OK: returncode normal, no FVS ERROR codes, summary rows "
            "consistent.  The guard is not simply rejecting everything.")
        say("   non-failing warnings reported: %s"
            % sorted({c["code"] for c in res["warnings"]}))
    except G.FVSRunError as exc:
        say("UNEXPECTED post-run failure: %s" % exc)
    proof["boundary40_passed"] = bool(pre40 and b40)
    proof["boundary40_summary_rows"] = r40["rows"]
    proof["boundary40_returncode"] = r40["rc"]
    proof["boundary40_codes"] = sorted({"%s:%s" % (c["code"], c["severity"])
                                        for c in codes40})
    proof["note_fvs03_warning"] = (
        "FVS03 is a WARNING (forest code outside model range, default used) "
        "emitted by the synthetic stand on both the 40 and 41 cases. "
        "check_run fails on ERROR severity only, so FVS03 does not fail a "
        "legal run.  This was found by the boundary test, not by inspection.")

    # ---- verdict -----------------------------------------------------------
    say("")
    say("=" * 78)
    say("VERDICT")
    say("-" * 78)
    say("MAXCYC parsed from PRGPRM.f90          : %d" % proof["maxcyc_parsed"])
    say("unguarded NUMCYCLE 41 return code      : %s" % proof["unguarded_returncode"])
    say("unguarded NUMCYCLE 41 summary rows     : %s" % proof["unguarded_summary_rows"])
    say("unguarded FVS04 listing line           : %s" % proof["unguarded_fvs04_line"])
    say("guard rejected 41 before launch        : %s" % proof["guarded_presubmit_rejected"])
    say("guard caught 41 after the run          : %s" % proof["guarded_postrun_detected"])
    say("boundary NUMCYCLE 40 passed guarded    : %s" % proof["boundary40_passed"])
    say("boundary NUMCYCLE 40 summary rows      : %s" % proof["boundary40_summary_rows"])
    say("=" * 78)

    with open(os.path.join(WORK, "PROOF.txt"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    with open(os.path.join(WORK, "proof.json"), "w") as fh:
        json.dump(proof, fh, indent=2, sort_keys=True)
    print("\nwrote PROOF.txt and proof.json to " + WORK)


if __name__ == "__main__":
    main()
