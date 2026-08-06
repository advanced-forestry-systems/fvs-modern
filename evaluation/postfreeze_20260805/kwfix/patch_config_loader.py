#!/usr/bin/env python3
"""Post-freeze patch for the BAMAX / BAIMULT / HTGMULT keyword writers.

Reads the LIVE config/config_loader.py, writes TWO files into the post-freeze
working directory, and does NOT modify the live file:

    kwfix/config_loader.py_PREV_20260805   byte copy of the live (pre-fix) file
    kwfix/config_loader.py                 the corrected file

The live file stays at md5 abc36de49ab8eb458d4e8916e55c466e, which is the
provenance tag recorded in FVS_MODEL_FREEZE_2026-08-05.md Table 1 for the file
the frozen run imported at runtime.
"""
import hashlib
import os
import shutil
import sys

SRC = "/users/PUOM0008/crsfaaron/fvs-modern/config/config_loader.py"
DST_DIR = "/fs/scratch/PUOM0008/crsfaaron/postfreeze_20260805/kwfix"
FROZEN_MD5 = "abc36de49ab8eb458d4e8916e55c466e"

OLD_BAMAX = '''    def _format_bamax_keywords(self, values: list, comments: bool) -> str:
        """Format BAMAX keyword block."""
        lines = []
        if comments:
            lines.append("!! Maximum basal area")
        # BAMAX uses a single value or per species
        if isinstance(values, (int, float)):
            lines.append(f"BAMAX           {values:10.1f}")
        else:
            for i, val in enumerate(values):
                if isinstance(val, str) or val is None:
                    continue
                if val > 0:
                    lines.append(f"BAMAX           {i + 1:10d}{val:10.1f}")
        return "\\n".join(lines)
'''

NEW_BAMAX = '''    def _format_bamax_keywords(self, values: list, comments: bool) -> str:
        """Format BAMAX keyword block.

        COLUMN LAYOUT, read out of src-converted/vbase/initre.f90 option 66
        ("OPTION NUMBER 66: BAMAX") and then proved by FVS keyword echo against
        the frozen binary on 2026-08-05:

            6800 CONTINUE
            LMORT=.TRUE.
            IF(ARRAY(1).GT.0.0) THEN
              BAMAX=ARRAY(1)
              LBAMAX=.TRUE.
            ENDIF

        BAMAX is NOT a scheduled activity.  It carries no date field and no
        species field; the maximum basal area is read straight out of field 1,
        which keyrdr.f90 takes from RECORD(11:20).  The previous writer padded
        the keyword to 16 columns, which pushed the value into columns 17-26 and
        left field 1 blank, so ARRAY(1) was 0.0, the .GT.0.0 test failed, and
        BAMAX was never set.  FVS raised no error, because a blank field is
        perfectly legal.  Proof: both legacy forms echoed
        "BAMAX      MAXIMUM BASAL AREA=      0.00" while the corrected form
        echoed 250.00, in postfreeze_20260805/kwfix/case_bamax/run.out.

        The previous per-species branch was doubly wrong.  BAMAX has no species
        dimension at all, so a per-species vector cannot be expressed in this
        keyword.  It is now rejected loudly instead of written as records FVS
        discards without comment.
        """
        lines = []
        if comments:
            lines.append("!! Maximum basal area (ft2 ac-1), stand level, no species field")
        if isinstance(values, bool):
            raise TypeError("BAMAX value must be numeric, got bool")
        if isinstance(values, (int, float)):
            if float(values) > 0:
                lines.append(f"{'BAMAX':<10}{float(values):10.1f}")
            return "\\n".join(lines)
        numeric = [
            float(v) for v in values
            if v is not None and not isinstance(v, str) and float(v) > 0
        ]
        distinct = sorted({round(v, 4) for v in numeric})
        if not distinct:
            return "\\n".join(lines)
        if len(distinct) > 1:
            raise ValueError(
                "BAMAX is a stand-level keyword with no species field "
                "(initre.f90 option 66 reads only ARRAY(1)), but "
                f"{len(distinct)} distinct per-species values were supplied "
                f"(first five: {distinct[:5]}). Collapse them to a single stand "
                "value, or use SDIMAX, which does take a species index."
            )
        lines.append(f"{'BAMAX':<10}{distinct[0]:10.1f}")
        return "\\n".join(lines)
'''

OLD_BAIMULT = '''    def _format_baimult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format BAIMULT (growth multiplier) keyword block."""
        lines = []
        if comments:
            lines.append("!! Diameter growth multipliers (calibrated / default)")
        for i, mult in enumerate(multipliers):
            if abs(mult - 1.0) > 0.01:
                # BAIMULT via READCORD or GROWTH multiplier approach
                # Using species level growth multiplier
                lines.append(f"BAIMULT         {i + 1:10d}{mult:10.4f}")
        return "\\n".join(lines)
'''

NEW_BAIMULT = '''    def _format_scheduled_multiplier(
        self, keyword: str, multipliers, comments: bool, header: str
    ) -> str:
        """Emit a scheduled per-species multiplier keyword block.

        Shared by BAIMULT (initre.f90 option 58, activity code 91) and HTGMULT
        (option 62, activity code 92).  HTGMULT is literally "6400 CONTINUE;
        I=92; GOTO 6005", i.e. it re-enters the BAIMULT handler, so the two have
        identical record layouts.  That handler reads:

            IDT=1
            IF (LNOTBK(1)) IDT=IFIX(ARRAY(1))   <- field 1 is the DATE
            IF (.NOT.LNOTBK(3)) ARRAY(3)=1.0
            CALL SPDECD (2,IS,...)              <- field 2 is the SPECIES
            CALL OPNEW(KODE,IDT,I,2,ARRAY(2))   <- field 3 is the MULTIPLIER

        Both are therefore SCHEDULED activities with the same layout as MORTMULT:
        cols 1-10 keyword, f1 cols 11-20 date (blank -> cycle 1), f2 cols 21-30
        species, f3 cols 31-40 multiplier.  The previous writer padded the
        keyword to 16 columns and omitted the date field, putting every value six
        columns right of where keyrdr.f90 reads it.

        That offset is VALUE DEPENDENT, which is exactly why it survived so long.
        Proved by FVS keyword echo against the frozen binary on 2026-08-05, in
        postfreeze_20260805/kwfix/case_mult/run.out: species 2 with multiplier
        1.8914 echoed identically under both layouts, because the six-column
        spill happened to land in blank columns, whereas species 12 with
        multiplier 12.3456 raised "FVS04 ERROR: A REQUIRED PARAMETER IS MISSING
        OR A PARAMETER IS INCORRECT; KEYWORD IGNORED" under the legacy layout and
        decoded correctly under this one.  The production driver discards stdout,
        stderr and the return code, so that FVS04 was invisible and the keyword
        was dropped in silence.

        Suppression is also announced.  A multiplier within 0.01 of unity is a
        no-op and is correctly omitted, but a family in which EVERY species is
        unity produces an empty block that is indistinguishable from a broken
        writer.  That ambiguity is what made the fvs_regional arm look
        mis-formatted when it was in fact being handed an all-unity vector.  When
        a non-empty vector yields zero records, a one-line note now goes to
        stderr, and to the block itself when comments are on.
        """
        lines = []
        if comments:
            lines.append(header)
        n_emitted = 0
        for i, mult in enumerate(multipliers):
            if abs(float(mult) - 1.0) > 0.01:
                lines.append(
                    f"{keyword:<10}{'':10}{i + 1:10d}{float(mult):10.4f}"
                )
                n_emitted += 1
        n_in = len(multipliers)
        if n_in and n_emitted == 0:
            msg = (
                f"[config_loader] {keyword}: {n_in} species supplied, 0 records "
                "emitted because every multiplier is within 0.01 of unity. This "
                "family is a NO-OP for this variant; the arm runs on the "
                "uncalibrated growth model. This is a calibration content "
                "question, not a keyword formatting failure."
            )
            sys.stderr.write(msg + "\\n")
            if comments:
                lines.append("!! " + msg[len("[config_loader] "):])
        return "\\n".join(lines)

    def _format_baimult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format BAIMULT (diameter growth multiplier) keyword block.

        See _format_scheduled_multiplier for the column layout and its proof.
        """
        return self._format_scheduled_multiplier(
            "BAIMULT", multipliers, comments,
            "!! Diameter growth multipliers (calibrated / default)",
        )
'''

OLD_HTGMULT = '''    def _format_htgmult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format height growth multiplier keyword block."""
        lines = []
        if comments:
            lines.append("!! Height growth multipliers (calibrated / default)")
        for i, mult in enumerate(multipliers):
            if abs(mult - 1.0) > 0.01:
                lines.append(f"HTGMULT         {i + 1:10d}{mult:10.4f}")
        return "\\n".join(lines)
'''

NEW_HTGMULT = '''    def _format_htgmult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format HTGMULT (height growth multiplier) keyword block.

        See _format_scheduled_multiplier for the column layout and its proof.
        """
        return self._format_scheduled_multiplier(
            "HTGMULT", multipliers, comments,
            "!! Height growth multipliers (calibrated / default)",
        )
'''


def main():
    src = open(SRC).read()
    md5 = hashlib.md5(src.encode()).hexdigest()
    print(f"live config_loader.py md5 = {md5}")
    if md5 != FROZEN_MD5:
        print(f"WARNING: live md5 differs from the freeze tag {FROZEN_MD5}")

    os.makedirs(DST_DIR, exist_ok=True)
    shutil.copy(SRC, os.path.join(DST_DIR, "config_loader.py_PREV_20260805"))

    out = src
    for name, old, new in (
        ("BAMAX", OLD_BAMAX, NEW_BAMAX),
        ("BAIMULT", OLD_BAIMULT, NEW_BAIMULT),
        ("HTGMULT", OLD_HTGMULT, NEW_HTGMULT),
    ):
        if old not in out:
            print(f"FAIL: exact block for {name} not found; aborting")
            sys.exit(1)
        if out.count(old) != 1:
            print(f"FAIL: block for {name} matched {out.count(old)} times; aborting")
            sys.exit(1)
        out = out.replace(old, new)
        print(f"patched {name}")

    # make sure `sys` is importable in the module.  It must go AFTER any
    # `from __future__` line, which Python requires to be the first statement.
    lines = out.split("\n")
    if not any(ln.strip() == "import sys" for ln in lines):
        anchor = None
        for idx, ln in enumerate(lines):
            if ln.startswith("from __future__ import"):
                anchor = idx + 1
        if anchor is None:
            for idx, ln in enumerate(lines):
                if ln.startswith("import ") or ln.startswith("from "):
                    anchor = idx
                    break
        lines.insert(anchor, "import sys")
        out = "\n".join(lines)
        print(f"inserted 'import sys' at line {anchor + 1}")
    else:
        print("'import sys' already present")

    dst = os.path.join(DST_DIR, "config_loader.py")
    with open(dst, "w") as fh:
        fh.write(out)
    print(f"wrote {dst}")
    print("new md5 =", hashlib.md5(out.encode()).hexdigest())
    print("live file untouched, md5 still", hashlib.md5(open(SRC, 'rb').read()).hexdigest())


if __name__ == "__main__":
    main()
