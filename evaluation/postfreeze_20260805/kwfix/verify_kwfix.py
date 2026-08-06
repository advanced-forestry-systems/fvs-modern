#!/usr/bin/env python3
"""Verify the post-freeze keyword-writer fix.

Three things are checked, all against real objects, none asserted:
 1. ZERO IMPACT ON THE FREEZE. The keyword block the driver actually generates
    for the NE calibrated config must be byte identical before and after the
    fix, for every arm's calib_components setting used in the frozen run.
 2. The corrected writers put values in the columns keyrdr.f90 reads.
 3. BAMAX with a genuine per-species vector is now rejected loudly.
"""
import hashlib
import importlib.util
import io
import json
import sys
import contextlib

import numpy as np

OLD = "/users/PUOM0008/crsfaaron/fvs-modern/config/config_loader.py"
NEW = "/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/kwfix/config_loader.py"
CONFIG_DIR = "/users/PUOM0008/crsfaaron/fvs-modern/config"
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-modern")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def block(loader_cls, include_dg, include_hg, include_mort):
    """Reproduce run_driver_final.partial_calibrated_keywords exactly."""
    ld = loader_cls("ne", version="calibrated", config_dir=CONFIG_DIR)
    cats = ld.config.get("categories", {})
    pre = ld.config.get("calibration_multipliers") or {}
    parts = []
    if bool(ld.config.get("_emit_sdimax", True)):
        sdi = ld._find_sdi_param(cats)
        if sdi is not None:
            parts.append(ld._format_sdimax_keywords(sdi, False))
    bam = ld._find_bamax_param(cats)
    if bam is not None:
        parts.append(ld._format_bamax_keywords(bam, False))
    if include_mort:
        mm = ld._array_or_none(pre.get("mort_multiplier"))
        if mm is None:
            mm = ld._compute_mortality_multipliers(cats)
        if mm is not None:
            parts.append(ld._format_mortmult_keywords(mm, False))
    if include_dg:
        dds = ld._array_or_none(pre.get("dds_multiplier"))
        gm = np.sqrt(np.clip(dds, 0.01, 100.0)) if dds is not None \
            else ld._compute_growth_multipliers(cats)
        if gm is not None:
            parts.append(ld._format_baimult_keywords(gm, False))
    if include_hg:
        hm = ld._array_or_none(pre.get("htg_multiplier"))
        if hm is None:
            hm = ld._compute_height_multipliers(cats)
        if hm is not None:
            parts.append(ld._format_htgmult_keywords(hm, False))
    return "\n".join(p for p in parts if p)


res = {}
log = []


def say(s=""):
    print(s)
    log.append(s)


mod_old = load(OLD, "cl_old")
mod_new = load(NEW, "cl_new")
COld, CNew = mod_old.FvsConfigLoader, mod_new.FvsConfigLoader

say("=== 1. Zero impact on the frozen run ===")
say("The frozen five-arm run used calib_components None (fvs_base), "
    "{dg,hg,mort all True} (fvs_regional), and all False (three Greg arms).")
same = True
for label, (dg, hg, mo) in {
    "fvs_regional (dg=T,hg=T,mort=T)": (True, True, True),
    "greg arms (dg=F,hg=F,mort=F)": (False, False, False),
}.items():
    err_o, err_n = io.StringIO(), io.StringIO()
    with contextlib.redirect_stderr(err_o):
        b_old = block(COld, dg, hg, mo)
    with contextlib.redirect_stderr(err_n):
        b_new = block(CNew, dg, hg, mo)
    h_old = hashlib.md5(b_old.encode()).hexdigest()
    h_new = hashlib.md5(b_new.encode()).hexdigest()
    ident = b_old == b_new
    same &= ident
    say(f"  {label}")
    say(f"    records old={len(b_old.splitlines())} new={len(b_new.splitlines())}")
    say(f"    md5 old={h_old} new={h_new}  IDENTICAL={ident}")
    if err_n.getvalue():
        say("    new-writer stderr note(s):")
        for ln in err_n.getvalue().strip().splitlines():
            say("      " + ln)
    res[f"identical::{label}"] = ident
res["freeze_zero_impact"] = bool(same)
say(f"  ZERO IMPACT ON FREEZE = {same}")
say()

say("=== 2. Corrected column layout ===")
ld = CNew("ne", version="calibrated", config_dir=CONFIG_DIR)


def fields(rec):
    r = rec.ljust(60)
    return {"kw": r[0:8], "f1": r[10:20], "f2": r[20:30], "f3": r[30:40]}


mult = np.ones(108)
mult[1] = 1.8914
mult[11] = 12.3456
for kw, fn in (("BAIMULT", ld._format_baimult_keywords),
               ("HTGMULT", ld._format_htgmult_keywords)):
    blk = fn(mult, False)
    say(f"  {kw}:")
    for rec in blk.splitlines():
        say(f"    |{rec}|  {json.dumps(fields(rec))}")
    res[f"{kw}_records"] = blk.splitlines()

say("  BAMAX scalar 250.0:")
rec = ld._format_bamax_keywords(250.0, False)
say(f"    |{rec}|  {json.dumps(fields(rec))}")
res["BAMAX_scalar"] = rec
say()

say("=== 3. BAMAX per-species vector is now rejected ===")
try:
    ld._format_bamax_keywords([200.0, 250.0, 300.0], False)
    say("  NOT REJECTED -- FAIL")
    res["bamax_vector_rejected"] = False
except ValueError as e:
    say(f"  ValueError: {e}")
    res["bamax_vector_rejected"] = True
say("  BAMAX uniform vector collapses to one record:")
say(f"    |{ld._format_bamax_keywords([250.0, 250.0, 250.0], False)}|")
say()

say("=== 4. Suppression diagnostic (defect 3 self-announcement) ===")
err = io.StringIO()
with contextlib.redirect_stderr(err):
    blk = ld._format_baimult_keywords(np.ones(108), False)
say(f"  emitted records = {len(blk.splitlines()) if blk.strip() else 0}")
say(f"  emitted bytes   = {len(blk)}  (must be 0 for driver-path safety)")
say("  stderr:")
for ln in err.getvalue().strip().splitlines():
    say("    " + ln)
res["suppression_note_emitted"] = bool(err.getvalue().strip())
res["suppression_block_bytes"] = len(blk)

with open("/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/kwfix/PROOF_kwfix.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
with open("/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/kwfix/proof_kwfix.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
print("\nwrote PROOF_kwfix.txt")
