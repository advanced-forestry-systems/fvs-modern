#!/usr/bin/env python3
"""
make_guarded_driver.py -- produce run_driver_guarded.py from run_driver_final.py.

The original is NEVER modified in place: its md5 (6072753c4c069fc273e11347d199a7a2)
is a published provenance tag.  This script reads it, applies three surgical
replacements, and writes the guarded COPY into the maxcyc working directory.
"""
import hashlib
import os
import sys

SRC = "/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804/run_driver_final.py"
WORK = "/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/maxcyc"
DST = os.path.join(WORK, "run_driver_guarded.py")
EXPECT_MD5 = "6072753c4c069fc273e11347d199a7a2"

src = open(SRC).read()
md5 = hashlib.md5(src.encode()).hexdigest()
print("source md5: %s (expected %s) match=%s"
      % (md5, EXPECT_MD5, md5 == EXPECT_MD5))
if md5 != EXPECT_MD5:
    sys.exit("REFUSING: source md5 does not match the published provenance tag")

# ---------------------------------------------------------------------------
# 1. Import the guard module and resolve MAXCYC from source.
# ---------------------------------------------------------------------------
A_OLD = 'FIA_DIR      = "/fs/scratch/PUOM0008/crsfaaron/FIA"\n'
A_NEW = '''FIA_DIR      = "/fs/scratch/PUOM0008/crsfaaron/FIA"

# --- MAXCYC / NUMCYCLE SILENT-TRUNCATION GUARD (2026-08-05) -----------------
# Diagnostic only: changes no model output, requires no rebuild.
# FVS is compiled with PARAMETER (MAXCYC=40) in src-converted/base/PRGPRM.f90.
# NUMCYCLE > MAXCYC does NOT abort.  FVS emits "FVS04 ERROR: A REQUIRED
# PARAMETER IS MISSING OR A PARAMETER IS INCORRECT; KEYWORD IGNORED" into the
# .out listing, discards the keyword, runs a SINGLE cycle, exits with return
# code 20, and still writes a well-formed TWO-ROW FVS_Summary2 table.  This
# driver used to send stdout and stderr to DEVNULL and never inspected the
# return code, so the failure was completely invisible in production.
sys.path.insert(0, os.path.join(PROJECT_ROOT, "tools"))
from fvs_run_guards import (  # noqa: E402
    FVSKeywordError, FVSRunError, check_run, read_maxcyc, validate_numcycle,
)

MAXCYC = read_maxcyc(os.path.join(PROJECT_ROOT, "src-converted", "base",
                                  "PRGPRM.f90"))
'''

# ---------------------------------------------------------------------------
# 2. Pre-submit guard + capture stdout/stderr instead of discarding them.
# ---------------------------------------------------------------------------
B_OLD = '''        key = os.path.join(tmp, "run.key")
        with open(key, "w") as fh:
            fh.write(KEYFILE.format(sid=sid, db=db,
                                    cycle_kw=build_cycle_kw(ncyc),
                                    calib_kw=calib_kw, htg_kw=htg_kw))

        subprocess.run(
            [binary, f"--keywordfile={key}"],
            cwd=tmp,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=env,
            timeout=120,
        )
'''
B_NEW = '''        key = os.path.join(tmp, "run.key")
        # PRE-SUBMIT GUARD: refuse to launch a run FVS would silently truncate.
        # In "annual" mode NUMCYCLE carries ncyc and is the value at risk; in
        # "single" mode NUMCYCLE is fixed at 1.
        n_cycles_req = ncyc if CYCLE_MODE == "annual" else 1
        validate_numcycle(n_cycles_req, maxcyc=MAXCYC)
        with open(key, "w") as fh:
            fh.write(KEYFILE.format(sid=sid, db=db,
                                    cycle_kw=build_cycle_kw(ncyc),
                                    calib_kw=calib_kw, htg_kw=htg_kw))

        # Capture stdout and stderr instead of discarding them to DEVNULL, and
        # keep the process object so the return code can be inspected.
        proc = subprocess.run(
            [binary, f"--keywordfile={key}"],
            cwd=tmp,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            timeout=120,
        )
'''

# ---------------------------------------------------------------------------
# 3. Post-run guard before returning, and loud failure on a guard trip.
# ---------------------------------------------------------------------------
C_OLD = '''        con.close()
        return df, tl, None

    except subprocess.TimeoutExpired:
        return None, None, "timeout"
'''
C_NEW = '''        con.close()

        # POST-RUN GUARD: return code, FVS error codes in the listing, and
        # summary-row / cycle-count consistency.  Raises FVSRunError rather
        # than handing back a plausible but truncated two-row summary.
        check_run(proc.returncode, os.path.join(tmp, "run.out"),
                  expected_cycles=n_cycles_req,
                  summary_rows=(0 if df is None else len(df)),
                  stderr_text=proc.stderr.decode("utf8", "replace"))
        return df, tl, None

    except (FVSKeywordError, FVSRunError) as exc:
        # Fail LOUDLY.  Do not degrade a guard trip into a soft error string
        # that a caller can ignore; SystemExit is not caught by the generic
        # "except Exception" below.
        print("FATAL FVS GUARD [%s variant=%s model=%s ncyc=%s]: %s"
              % (sid, variant, model_name, ncyc, exc), file=sys.stderr)
        raise SystemExit(3)
    except subprocess.TimeoutExpired:
        return None, None, "timeout"
'''

out = src
for i, (old, new) in enumerate([(A_OLD, A_NEW), (B_OLD, B_NEW), (C_OLD, C_NEW)], 1):
    n = out.count(old)
    if n != 1:
        sys.exit("REFUSING: replacement %d matched %d times, expected 1" % (i, n))
    out = out.replace(old, new)
    print("applied replacement %d" % i)

with open(DST, "w") as fh:
    fh.write(out)
print("wrote %s" % DST)
