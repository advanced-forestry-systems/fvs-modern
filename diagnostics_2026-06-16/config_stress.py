#!/usr/bin/env python3
"""
Full-repo config stress test (2026-06-16). For every FVS variant with a calibrated
config and executable, load the config, emit generate_keywords (post-fix loader),
and tally the keyword blocks. Confirms the configs are valid and the SDIMAX
emission is suppressed (0 lines) after the WO-1 fix, while the multiplier blocks
remain. Pure load/emit; no projection, no region mapping needed.
"""
import os, sys, glob, json, traceback

P = "/users/PUOM0008/crsfaaron/fvs-modern"
WORK = os.path.expanduser("~/overthin_work")
os.environ["FVS_PROJECT_ROOT"] = P
os.environ["FVS_CONFIG_DIR"] = P + "/config"
for p in [WORK, P, P + "/calibration/python"]:
    sys.path.insert(0, p)
from config.config_loader import FvsConfigLoader

CFG = P + "/config"
LIB = P + "/lib"


def nonunity(a):
    return sum(1 for x in (a or []) if isinstance(x, (int, float)) and abs(x - 1.0) > 1e-9)


variants = sorted(os.path.basename(f)[:-5] for f in glob.glob(CFG + "/calibrated/[a-z][a-z].json"))
print(f"{'var':4s} {'exe':3s} {'load':4s} {'SDIMAX':6s} {'MORT':5s} {'BAI':4s} {'HTG':4s} "
      f"{'dds_nu':6s} {'mort_nu':7s} {'htg_nu':6s} note")
n_ok = 0
n_fail = 0
for v in variants:
    has_exe = os.path.exists(os.path.join(LIB, "FVS" + v))
    try:
        ld = FvsConfigLoader(v, version="calibrated", config_dir=CFG)
        kw = ld.generate_keywords(include_comments=False)
        cnt = {k: sum(1 for l in kw.splitlines() if l.strip().startswith(k))
               for k in ("SDIMAX", "MORTMULT", "BAIMULT", "HTGMULT")}
        cm = ld.config.get("calibration_multipliers", {}) or {}
        dds_nu = nonunity(cm.get("dds_multiplier"))
        mort_nu = nonunity(cm.get("mort_multiplier"))
        htg_nu = nonunity(cm.get("htg_multiplier"))
        note = "OK" if cnt["SDIMAX"] == 0 else "!! SDIMAX still emitted"
        print(f"{v:4s} {'Y' if has_exe else '-':3s} {'OK':4s} "
              f"{cnt['SDIMAX']:<6d} {cnt['MORTMULT']:<5d} {cnt['BAIMULT']:<4d} {cnt['HTGMULT']:<4d} "
              f"{dds_nu:<6d} {mort_nu:<7d} {htg_nu:<6d} {note}")
        n_ok += 1
    except Exception as e:
        print(f"{v:4s} {'Y' if has_exe else '-':3s} FAIL  -- load/emit error: {e!r}")
        n_fail += 1
print(f"\nSUMMARY: {n_ok} configs load+emit OK, {n_fail} failed, of {len(variants)} calibrated configs.")
print("DONE_STRESS")
