#!/usr/bin/env python3
"""
fvs_run_guards.py -- pre-submit and post-run guards against FVS silent failure.

Created 2026-08-05.  Diagnostic only: this module changes no model output and
requires no rebuild.  It exists to make an already-existing failure visible.

THE DEFECT THIS GUARDS
----------------------
FVS is compiled with a fixed cycle-count ceiling:

    src-converted/base/PRGPRM.f90
        PARAMETER (MAXCYC=40)

A keyword file that asks for more than MAXCYC cycles, e.g.

    TIMEINT            0         1
    NUMCYCLE          41

does NOT abort.  It parses, the keyword handler emits

    FVS04 ERROR: A REQUIRED PARAMETER IS MISSING OR A PARAMETER IS INCORRECT;
    KEYWORD IGNORED

into the .out listing, the NUMCYCLE keyword is then discarded, FVS falls back to
a SINGLE cycle, exits with return code 20 (the normal FVS stop is 10), and
writes a perfectly well-formed TWO-ROW FVS_Summary2 table.  Every downstream
consumer sees valid-looking output.

The production driver (run_driver_final.py, run_stand()) sends both stdout and
stderr to subprocess.DEVNULL and never inspects the return code, so in
production this failure is completely invisible.  It cost 45 runs in the
2026-08-05 stress battery before the boundary at NUMCYCLE 40/41 was located.

WHAT THIS MODULE PROVIDES
-------------------------
    read_maxcyc(prgprm_path)        parse MAXCYC out of PRGPRM.f90 (not hardcoded)
    validate_numcycle(n, maxcyc)    PRE-SUBMIT guard  -> FVSKeywordError
    validate_keyword_file(path, ..) PRE-SUBMIT guard on a whole .key file
    scan_listing(out_path)          find FVS\\d\\d error codes in the .out listing
    check_run(rc, out_path, ...)    POST-RUN guard    -> FVSRunError

Standard library only.  pandas is optional and import is guarded, because the
driver environment may be minimal.

CLI:
    python3 fvs_run_guards.py --key run.key --out run.out --rc 20
"""
from __future__ import annotations

import argparse
import os
import re
import sys

__all__ = [
    "MAXCYC_DEFAULT",
    "FVS_NORMAL_RETURNCODES",
    "FVSGuardError",
    "FVSKeywordError",
    "FVSRunError",
    "read_maxcyc",
    "validate_numcycle",
    "validate_keyword_file",
    "scan_listing",
    "check_run",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Fallback used only when PRGPRM.f90 cannot be read or parsed.  The authoritative
#: value is whatever read_maxcyc() finds in the source.
MAXCYC_DEFAULT = 40

#: Default location of PRGPRM.f90 inside the fvs-modern tree.  tools/ sits one
#: level below the repo root, so ../src-converted/base/PRGPRM.f90.
DEFAULT_PRGPRM_PATH = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "src-converted", "base", "PRGPRM.f90")
)

#: FVS exits 10 on a normal STOP.  0 is accepted too for robustness across
#: builds; anything else (notably 20) is a failure.
FVS_NORMAL_RETURNCODES = (0, 10)

_MAXCYC_RE = re.compile(
    r"^\s*PARAMETER\s*\(\s*MAXCYC\s*=\s*(\d+)\s*\)", re.IGNORECASE | re.MULTILINE
)

#: Any FVS diagnostic code of the form FVSnn (FVS04, FVS07, ...).  Deliberately
#: general so codes not yet seen are still caught.
_FVS_CODE_RE = re.compile(r"\bFVS(\d\d)\b")

_NUMCYCLE_RE = re.compile(r"^\s*NUMCYCLE\b(.*)$", re.IGNORECASE)
_TIMEINT_RE = re.compile(r"^\s*TIMEINT\b(.*)$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------
class FVSGuardError(Exception):
    """Base class for every guard failure raised by this module."""


class FVSKeywordError(FVSGuardError):
    """A keyword value would be silently rejected by FVS.  Raised BEFORE launch."""


class FVSRunError(FVSGuardError):
    """An FVS invocation failed or emitted an error code.  Raised AFTER the run."""


# ---------------------------------------------------------------------------
# MAXCYC discovery
# ---------------------------------------------------------------------------
def read_maxcyc(prgprm_path: str | None = None) -> int:
    """Parse MAXCYC out of base/PRGPRM.f90.

    Returns the integer compiled into the binary.  Falls back to MAXCYC_DEFAULT
    if the file is unreadable or the PARAMETER statement is not found, so the
    guard degrades to "still guards, just with the documented default" rather
    than crashing the driver.
    """
    path = prgprm_path or DEFAULT_PRGPRM_PATH
    try:
        with open(path, "r", errors="replace") as fh:
            src = fh.read()
    except OSError:
        return MAXCYC_DEFAULT
    m = _MAXCYC_RE.search(src)
    if not m:
        return MAXCYC_DEFAULT
    return int(m.group(1))


# ---------------------------------------------------------------------------
# PRE-SUBMIT guards
# ---------------------------------------------------------------------------
def validate_numcycle(numcycle, maxcyc: int | None = None) -> int:
    """PRE-SUBMIT guard on the NUMCYCLE value.

    Raises FVSKeywordError if numcycle is not an integer >= 1 and <= MAXCYC.
    Returns the validated integer so it can be used inline.
    """
    if maxcyc is None:
        maxcyc = read_maxcyc()
    try:
        n = int(numcycle)
    except (TypeError, ValueError):
        raise FVSKeywordError(
            "NUMCYCLE is not an integer: %r" % (numcycle,)
        ) from None
    if n < 1:
        raise FVSKeywordError(
            "NUMCYCLE=%d is below the legal minimum of 1; FVS would reject the "
            "keyword and silently run a single cycle." % n
        )
    if n > maxcyc:
        raise FVSKeywordError(
            "NUMCYCLE=%d exceeds the compiled ceiling MAXCYC=%d "
            "(PARAMETER (MAXCYC=%d) in src-converted/base/PRGPRM.f90). "
            "FVS would NOT abort: it emits 'FVS04 ERROR: A REQUIRED PARAMETER "
            "IS MISSING OR A PARAMETER IS INCORRECT; KEYWORD IGNORED' into the "
            ".out listing, discards the keyword, runs ONE cycle, exits with "
            "return code 20, and still writes a well-formed two-row "
            "FVS_Summary2 table.  Refusing to launch."
            % (n, maxcyc, maxcyc)
        )
    return n


def _first_int(tail: str):
    """First integer token on the remainder of a fixed-field keyword line."""
    m = re.search(r"-?\d+", tail)
    return int(m.group(0)) if m else None


def validate_keyword_file(path: str, maxcyc: int | None = None) -> dict:
    """PRE-SUBMIT guard on a whole .key file, applied before FVS is invoked.

    Scans for NUMCYCLE (and TIMEINT if present), validates NUMCYCLE against
    MAXCYC, and returns a dict of what it found:
        {"numcycle": int|None, "timeint": int|None, "maxcyc": int, "path": str}
    Raises FVSKeywordError on an illegal value.
    """
    if maxcyc is None:
        maxcyc = read_maxcyc()
    if not os.path.exists(path):
        raise FVSKeywordError("keyword file does not exist: %s" % path)

    numcycle = None
    timeint = None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            m = _NUMCYCLE_RE.match(line)
            if m:
                v = _first_int(m.group(1))
                if v is not None:
                    numcycle = v
                continue
            m = _TIMEINT_RE.match(line)
            if m:
                # TIMEINT field 1 is the cycle index (0 = all), field 2 is the
                # cycle length in years.  The length is the value of interest.
                ints = [int(x) for x in re.findall(r"-?\d+", m.group(1))]
                if len(ints) >= 2:
                    timeint = ints[1]
                elif ints:
                    timeint = ints[0]

    if numcycle is not None:
        validate_numcycle(numcycle, maxcyc=maxcyc)
    if timeint is not None and timeint < 0:
        raise FVSKeywordError(
            "TIMEINT cycle length %d is negative in %s" % (timeint, path)
        )
    return {"numcycle": numcycle, "timeint": timeint,
            "maxcyc": maxcyc, "path": path}


# ---------------------------------------------------------------------------
# Listing scan
# ---------------------------------------------------------------------------
def _severity(line: str) -> str:
    """Classify an FVSnn diagnostic line as ERROR, WARNING, or OTHER.

    FVS tags its listing diagnostics with a severity word on the same line,
    e.g. 'FVS04 ERROR: ...' versus 'FVS03 WARNING: ...'.  The distinction
    matters: a legal NUMCYCLE 40 run on a synthetic stand still emits FVS03
    WARNING (forest code outside the model range, default used), which is
    informational and must NOT fail the run.  Only ERROR is a hard stop.
    """
    u = line.upper()
    if "ERROR" in u:
        return "ERROR"
    if "WARNING" in u or "CAUTION" in u:
        return "WARNING"
    return "OTHER"


def scan_listing(out_path: str) -> list:
    """Read the FVS .out listing and return every FVSnn diagnostic found.

    Returns a list of dicts:
        {"code": "FVS04", "severity": "ERROR", "lineno": int, "line": str}

    The FVS\\d\\d pattern is deliberately general so codes not yet encountered
    are still caught; severity is parsed from the same line so callers can
    separate hard errors from informational warnings.

    An unreadable or absent listing returns [] rather than raising.
    """
    results = []
    try:
        with open(out_path, "r", errors="replace") as fh:
            for i, line in enumerate(fh, start=1):
                seen = set()
                for m in _FVS_CODE_RE.finditer(line):
                    code = "FVS" + m.group(1)
                    if code in seen:
                        continue
                    seen.add(code)
                    results.append({
                        "code": code,
                        "severity": _severity(line),
                        "lineno": i,
                        "line": line.rstrip("\n").rstrip(),
                    })
    except OSError:
        return []
    return results


# ---------------------------------------------------------------------------
# POST-RUN guard
# ---------------------------------------------------------------------------
def check_run(returncode,
              out_path: str | None = None,
              expected_cycles: int | None = None,
              summary_rows: int | None = None,
              *,
              stderr_text: str | None = None,
              allow_codes=(),
              fail_on_severity=("ERROR", "OTHER")) -> dict:
    """POST-RUN guard.  Raises FVSRunError on any sign of the silent failure.

    Fails if:
      * returncode is not in FVS_NORMAL_RETURNCODES (10 = normal FVS STOP);
      * an FVSnn diagnostic of failing severity appears in the .out listing.
        By default ERROR and unclassified (OTHER) codes fail; WARNING codes are
        reported but do not fail, because a perfectly legal run emits FVS03
        WARNING on a synthetic stand.  Widen or narrow with fail_on_severity,
        or exempt specific codes with allow_codes;
      * summary_rows is inconsistent with expected_cycles.  FVS writes one
        summary row per cycle boundary, so N cycles yields N+1 rows.  The
        signature of this defect is expected_cycles=41 but summary_rows=2.

    Returns a dict of what it checked when everything passes, including any
    non-failing warnings so the caller can still log them.
    """
    problems = []
    codes = scan_listing(out_path) if out_path else []
    allow = set(allow_codes)
    fail_sev = set(fail_on_severity)
    warnings = [c for c in codes
                if c["code"] in allow or c["severity"] not in fail_sev]

    if returncode is not None and int(returncode) not in FVS_NORMAL_RETURNCODES:
        problems.append(
            "FVS returned code %s (normal FVS STOP is 10; 20 is the "
            "keyword-rejected / truncated-run signature)." % returncode
        )

    bad = [c for c in codes
           if c["code"] not in allow and c["severity"] in fail_sev]
    if bad:
        shown = "; ".join(
            "%s %s at listing line %d: %s"
            % (c["code"], c["severity"], c["lineno"], c["line"].strip())
            for c in bad[:5]
        )
        problems.append(
            "FVS error code(s) present in the listing %s -- %s"
            % (out_path, shown)
        )

    if expected_cycles is not None and summary_rows is not None:
        want = int(expected_cycles) + 1
        if int(summary_rows) != want:
            problems.append(
                "FVS_Summary2 has %d row(s) but %d cycle(s) were requested, so "
                "%d row(s) were expected.  A truncated run typically yields "
                "exactly 2 rows." % (summary_rows, expected_cycles, want)
            )

    if problems:
        msg = ("FVS run guard FAILED (returncode=%s, listing=%s):\n  - %s"
               % (returncode, out_path, "\n  - ".join(problems)))
        if stderr_text:
            msg += "\n  stderr: %s" % stderr_text.strip()[:500]
        raise FVSRunError(msg)

    return {"returncode": returncode, "out_path": out_path,
            "fvs_codes": [], "warnings": warnings,
            "expected_cycles": expected_cycles,
            "summary_rows": summary_rows, "ok": True}


# ---------------------------------------------------------------------------
# Optional pandas helper (guarded import; the driver env may be minimal)
# ---------------------------------------------------------------------------
def summary_row_count(db_path: str, table: str = "FVS_Summary2"):
    """Row count of the summary table, or None if it cannot be read.

    Uses sqlite3 from the standard library.  pandas is used only if available
    and is never required.
    """
    try:
        import sqlite3
    except ImportError:  # pragma: no cover
        return None
    try:
        con = sqlite3.connect(db_path)
        try:
            n = con.execute("SELECT COUNT(*) FROM %s" % table).fetchone()[0]
        finally:
            con.close()
        return int(n)
    except Exception:
        pass
    try:  # optional pandas path, guarded
        import pandas as pd  # noqa: F401
        import sqlite3
        con = sqlite3.connect(db_path)
        try:
            return int(len(pd.read_sql_query("SELECT * FROM %s" % table, con)))
        finally:
            con.close()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="fvs_run_guards.py",
        description="Pre-submit and post-run guards against FVS MAXCYC "
                    "silent truncation.")
    ap.add_argument("--key", help="FVS keyword (.key) file to validate")
    ap.add_argument("--out", help="FVS listing (.out) file to scan")
    ap.add_argument("--rc", type=int, default=None,
                    help="FVS process return code to check")
    ap.add_argument("--numcycle", type=int, default=None,
                    help="NUMCYCLE value to validate directly")
    ap.add_argument("--expected-cycles", type=int, default=None)
    ap.add_argument("--summary-rows", type=int, default=None)
    ap.add_argument("--db", default=None,
                    help="FVS_Data.db to count FVS_Summary2 rows from")
    ap.add_argument("--prgprm", default=None,
                    help="path to base/PRGPRM.f90 (default: repo-relative)")
    args = ap.parse_args(argv)

    maxcyc = read_maxcyc(args.prgprm)
    print("MAXCYC parsed from %s: %d"
          % (args.prgprm or DEFAULT_PRGPRM_PATH, maxcyc))

    status = 0

    if args.numcycle is not None:
        try:
            validate_numcycle(args.numcycle, maxcyc=maxcyc)
            print("PRE-SUBMIT OK: NUMCYCLE=%d <= MAXCYC=%d"
                  % (args.numcycle, maxcyc))
        except FVSKeywordError as exc:
            print("PRE-SUBMIT REJECT: %s" % exc)
            status = 2

    if args.key:
        try:
            info = validate_keyword_file(args.key, maxcyc=maxcyc)
            print("PRE-SUBMIT OK: %s (NUMCYCLE=%s TIMEINT=%s)"
                  % (args.key, info["numcycle"], info["timeint"]))
        except FVSKeywordError as exc:
            print("PRE-SUBMIT REJECT: %s" % exc)
            status = 2

    rows = args.summary_rows
    if rows is None and args.db:
        rows = summary_row_count(args.db)
        if rows is not None:
            print("FVS_Summary2 rows in %s: %d" % (args.db, rows))

    if args.out or args.rc is not None:
        if args.out:
            for c in scan_listing(args.out):
                print("LISTING %s %s line %d: %s"
                      % (c["code"], c["severity"], c["lineno"],
                         c["line"].strip()))
        try:
            check_run(args.rc, args.out,
                      expected_cycles=args.expected_cycles,
                      summary_rows=rows)
            print("POST-RUN OK: returncode=%s, no FVS error codes, "
                  "summary rows consistent." % args.rc)
        except FVSRunError as exc:
            print("POST-RUN FAIL: %s" % exc)
            status = 3

    return status


if __name__ == "__main__":
    sys.exit(_main())
