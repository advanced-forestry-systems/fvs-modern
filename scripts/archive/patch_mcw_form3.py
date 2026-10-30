#!/usr/bin/env python3
"""
Extend the Greg HG hook's MCW mechanism from a 4-slot GMCW to a 6-slot GMCW
carrying two per-species caps:
  slot 5 = CWMAX  (feet, 0 = no cap)  -> clamp MCW after evaluation
  slot 6 = DBHMAX (inches, 0 = no cap) -> evaluate the equation at
           DBHeff = min(DBH, DBHMAX) so quadratics do not turn over and
           linears do not run away.
New CSV schema: SPCD,FORM,A,B,C,CWMAX,DBHMAX (greg_mcw_coefficients.csv).
FORM 1: mcw = A + B*dbh + C*dbh^2 ; FORM 2: mcw = A * dbh^B. The >=0 floor and
FORM 1/2 branch are preserved.

Idempotent, anchored, need()/count guarded. Run ON CARDINAL from the fvs-modern
root (the dir that contains src-converted/). NOTE: the MCW loader lives in
base/greghghg.f90 (the HG hook), NOT base/gregdghg.f90 (that is the DG loader).
"""
import io, sys, os

def rw(path, fn):
    with io.open(path, "r", newline="") as f:
        s = f.read()
    out = fn(s)
    if out == s:
        print("  unchanged (already patched):", path)
        return
    with io.open(path, "w", newline="") as f:
        f.write(out)
    print("  patched:", path)

def need(s, tok, path, n=1):
    if s.count(tok) != n:
        sys.exit("ABORT: token count %d (want %d) in %s:\n  %r" % (s.count(tok), n, path, tok))

GREGMC = "src-converted/common/GREGMC.f90"
def p_gregmc(s):
    if "GMCW(MAXSP,6)" in s:
        return s
    need(s, "GMCW(MAXSP,4)", GREGMC)
    need(s, "!    GMCW      -- GMCW(ISPC,1:4): (1)=FORM code, (2:4)=A,B,C for MCW quadratic.", GREGMC)
    s = s.replace("GMCW(MAXSP,4)", "GMCW(MAXSP,6)", 1)
    s = s.replace(
        "!    GMCW      -- GMCW(ISPC,1:4): (1)=FORM code, (2:4)=A,B,C for MCW quadratic.",
        "!    GMCW      -- GMCW(ISPC,1:6): (1)=FORM,(2:4)=A,B,C,(5)=CWMAX,(6)=DBHMAX.", 1)
    return s

# MCW loader is in the HG hook file greghghg.f90.
GREGHGHG = "src-converted/base/greghghg.f90"
def p_greghghg(s):
    if "GMCW(ISPC,6)=0.0" in s:
        return s
    need(s, "GMCW(ISPC,1)=0.0; GMCW(ISPC,2)=0.0; GMCW(ISPC,3)=0.0; GMCW(ISPC,4)=0.0", GREGHGHG)
    need(s, "      READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2\n", GREGHGHG)
    need(s, "      TB(NG,1)=REAL(NN); TB(NG,2)=C0; TB(NG,3)=C1; TB(NG,4)=C2\n", GREGHGHG)
    need(s, "            GMCW(ISPC,3)=TB(J,3); GMCW(ISPC,4)=TB(J,4)\n", GREGHGHG)
    s = s.replace(
        "GMCW(ISPC,1)=0.0; GMCW(ISPC,2)=0.0; GMCW(ISPC,3)=0.0; GMCW(ISPC,4)=0.0",
        "GMCW(ISPC,1)=0.0; GMCW(ISPC,2)=0.0; GMCW(ISPC,3)=0.0; GMCW(ISPC,4)=0.0; GMCW(ISPC,5)=0.0; GMCW(ISPC,6)=0.0", 1)
    s = s.replace(
        "      READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2\n",
        "      READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2, C3, C4\n", 1)
    s = s.replace(
        "      TB(NG,1)=REAL(NN); TB(NG,2)=C0; TB(NG,3)=C1; TB(NG,4)=C2\n",
        "      TB(NG,1)=REAL(NN); TB(NG,2)=C0; TB(NG,3)=C1; TB(NG,4)=C2; TB(NG,5)=C3; TB(NG,6)=C4\n", 1)
    s = s.replace(
        "            GMCW(ISPC,3)=TB(J,3); GMCW(ISPC,4)=TB(J,4)\n",
        "            GMCW(ISPC,3)=TB(J,3); GMCW(ISPC,4)=TB(J,4)\n"
        "            GMCW(ISPC,5)=TB(J,5); GMCW(ISPC,6)=TB(J,6)\n", 1)
    return s

HTGF = "src-converted/ne/htgf.f90"
def p_htgf(s):
    if "DBHEFF" in s:
        return s
    need(s, "REAL OGCA(MAXTRE), MCWX, MCWA, MCWB, MCWC\n", HTGF)
    need(s, "MFORM=NINT(GMCW(ISP(JJ),1)); MCWA=GMCW(ISP(JJ),2); MCWB=GMCW(ISP(JJ),3); MCWC=GMCW(ISP(JJ),4)", HTGF)
    need(s, "MFORM=1; MCWA=8.0; MCWB=1.5; MCWC=0.0", HTGF)
    need(s,
         "    IF (MFORM.EQ.2) THEN\n"
         "      MCWX = MCWA * DBH(JJ)**MCWB\n"
         "    ELSE\n"
         "      MCWX = MCWA + MCWB*DBH(JJ) + MCWC*DBH(JJ)*DBH(JJ)\n"
         "    ENDIF\n"
         "    IF (MCWX.LT.0.0) MCWX = 0.0\n", HTGF)
    s = s.replace(
        "REAL OGCA(MAXTRE), MCWX, MCWA, MCWB, MCWC\n",
        "REAL OGCA(MAXTRE), MCWX, MCWA, MCWB, MCWC\nREAL MCWMAX, MCWDMX, DBHEFF\n", 1)
    # NOTE: free-form (-std=legacy) enforces a 132-column limit, so the two new
    # assignments go on a fresh continuation line rather than appending inline.
    s = s.replace(
        "MFORM=NINT(GMCW(ISP(JJ),1)); MCWA=GMCW(ISP(JJ),2); MCWB=GMCW(ISP(JJ),3); MCWC=GMCW(ISP(JJ),4)",
        "MFORM=NINT(GMCW(ISP(JJ),1)); MCWA=GMCW(ISP(JJ),2); MCWB=GMCW(ISP(JJ),3); MCWC=GMCW(ISP(JJ),4)\n"
        "      MCWMAX=GMCW(ISP(JJ),5); MCWDMX=GMCW(ISP(JJ),6)", 1)
    s = s.replace(
        "MFORM=1; MCWA=8.0; MCWB=1.5; MCWC=0.0",
        "MFORM=1; MCWA=8.0; MCWB=1.5; MCWC=0.0; MCWMAX=0.0; MCWDMX=0.0", 1)
    s = s.replace(
        "    IF (MFORM.EQ.2) THEN\n"
        "      MCWX = MCWA * DBH(JJ)**MCWB\n"
        "    ELSE\n"
        "      MCWX = MCWA + MCWB*DBH(JJ) + MCWC*DBH(JJ)*DBH(JJ)\n"
        "    ENDIF\n"
        "    IF (MCWX.LT.0.0) MCWX = 0.0\n",
        "    DBHEFF = DBH(JJ)\n"
        "    IF (MCWDMX.GT.0.0 .AND. DBH(JJ).GT.MCWDMX) DBHEFF = MCWDMX\n"
        "    IF (MFORM.EQ.2) THEN\n"
        "      MCWX = MCWA * DBHEFF**MCWB\n"
        "    ELSE\n"
        "      MCWX = MCWA + MCWB*DBHEFF + MCWC*DBHEFF*DBHEFF\n"
        "    ENDIF\n"
        "    IF (MCWX.LT.0.0) MCWX = 0.0\n"
        "    IF (MCWMAX.GT.0.0 .AND. MCWX.GT.MCWMAX) MCWX = MCWMAX\n", 1)
    return s

def main():
    for p in (GREGMC, GREGHGHG, HTGF):
        if not os.path.exists(p):
            sys.exit("ABORT: missing " + p)
    rw(GREGMC, p_gregmc)
    rw(GREGHGHG, p_greghghg)
    rw(HTGF, p_htgf)
    print("DONE.")

if __name__ == "__main__":
    main()
