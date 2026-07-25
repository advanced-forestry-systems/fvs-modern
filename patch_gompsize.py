#!/usr/bin/env python3
"""
Add a diameter (size) term to Greg Johnson's GOMPIT survival model.

Deployed (size-blind) hazard:
    eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4
Experimental (size-aware) hazard, matching greg_mortality_coefficients_sizeaware.csv
(SPCD,n,b0,b1,b2,b3,b4,b5,...  where b5 is the ln-DBH slope, DBH in inches):
    eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4 + b5*ln(DBH)

Four changed regions:
  1. common/GOMPMC.f90     GB(MAXSP,5) -> GB(MAXSP,6)  (COMMON /GOMPMR/ grows one slot)
  2. base/gompmort.f90     GOMPLOAD  -> read a 6th coefficient (b5) per species,
                           line-buffered so 5-coef (old) files still parse and
                           GB(ISPC,6) falls back to 0.0.
  3. base/gompmort.f90     GOMPSURV  -> accept DBHV and add B5*LOG(MAX(DBHV,0.1)).
  4. vls/morts.f90         pass the tree diameter D=DBH(I) into the GOMPSURV call.

NOTE ON THE CALL SITE: the NE binary (FVSne) links ../vls/morts.f90 for its
mortality driver (see bin/FVSne_sourceList.txt); there is no ne/morts.f90. So the
one NE-relevant GOMPSURV call site is vls/morts.f90, and that is what we patch.

Idempotent, anchored, need()/count guarded. Run ON CARDINAL from the fvs-modern
root (the dir that contains src-converted/).
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

# ---------------------------------------------------------------------------
# 1. common/GOMPMC.f90 : GB(MAXSP,5) -> GB(MAXSP,6)
# ---------------------------------------------------------------------------
GOMPMC = "src-converted/common/GOMPMC.f90"
def p_gompmc(s):
    if "GB(MAXSP,6)" in s:
        return s
    need(s, "REAL GB(MAXSP,5), CCHT(MAXTRE)", GOMPMC)
    need(s, "!    GB     -- per-FVS-species coefficient table, GB(ISPC,1:5) = b0,b1,b2,b3,b4", GOMPMC)
    s = s.replace("REAL GB(MAXSP,5), CCHT(MAXTRE)",
                  "REAL GB(MAXSP,6), CCHT(MAXTRE)", 1)
    s = s.replace(
        "!    GB     -- per-FVS-species coefficient table, GB(ISPC,1:5) = b0,b1,b2,b3,b4",
        "!    GB     -- per-FVS-species coefficient table, GB(ISPC,1:6) = b0,b1,b2,b3,b4,b5\n"
        "!              (b5 = ln-DBH size slope; 0.0 for size-blind coefficient files).",
        1)
    return s

# ---------------------------------------------------------------------------
# 2 & 3. base/gompmort.f90 : GOMPLOAD reader + GOMPSURV size term
# ---------------------------------------------------------------------------
GOMPMORT = "src-converted/base/gompmort.f90"
def p_gompmort(s):
    already = ("B5*LOG(MAX(DBHV,0.1))" in s)
    if already:
        return s

    # -- GOMPLOAD: local coefficient table + scalars grow one slot --
    need(s, "REAL    TB(MXG,5)", GOMPMORT)
    need(s, "REAL B0,B1,B2,B3,B4\n", GOMPMORT)
    need(s, "INTEGER I, J, NG, IOS, U, IFIA, NN, ISPC\n", GOMPMORT)
    s = s.replace("REAL    TB(MXG,5)", "REAL    TB(MXG,6)", 1)
    s = s.replace("REAL B0,B1,B2,B3,B4\n", "REAL B0,B1,B2,B3,B4,B5\n", 1)
    s = s.replace("INTEGER I, J, NG, IOS, U, IFIA, NN, ISPC\n",
                  "INTEGER I, J, NG, IOS, IOS2, U, IFIA, NN, ISPC\n", 1)

    # -- GOMPLOAD: zero-init loop 1..5 -> 1..6 --
    need(s, "  DO J=1,5\n    GB(ISPC,J) = 0.0\n  ENDDO\n", GOMPMORT)
    s = s.replace("  DO J=1,5\n    GB(ISPC,J) = 0.0\n  ENDDO\n",
                  "  DO J=1,6\n    GB(ISPC,J) = 0.0\n  ENDDO\n", 1)

    # -- GOMPLOAD: line-buffered read of the (up to) 6th coefficient --
    old_read = (
        "  READ(U,*,IOSTAT=IOS) IFIA, NN, B0, B1, B2, B3, B4\n"
        "  IF (IOS.NE.0) GO TO 20\n"
        "  IF (NG.GE.MXG) GO TO 20\n"
        "  NG = NG + 1\n"
        "  GSPCD(NG) = IFIA\n"
        "  TB(NG,1)=B0; TB(NG,2)=B1; TB(NG,3)=B2; TB(NG,4)=B3; TB(NG,5)=B4\n"
    )
    new_read = (
        "  READ(U,'(A)',IOSTAT=IOS) LINE\n"
        "  IF (IOS.NE.0) GO TO 20\n"
        "  IF (LINE.EQ.' ') GO TO 10\n"
        "  ! Try 6 coefficients (size-aware b5); fall back to 5 (size-blind) with\n"
        "  ! b5=0. Internal (line) read never over-consumes the next species row.\n"
        "  B5 = 0.0\n"
        "  READ(LINE,*,IOSTAT=IOS2) IFIA, NN, B0, B1, B2, B3, B4, B5\n"
        "  IF (IOS2.NE.0) THEN\n"
        "    B5 = 0.0\n"
        "    READ(LINE,*,IOSTAT=IOS2) IFIA, NN, B0, B1, B2, B3, B4\n"
        "    IF (IOS2.NE.0) GO TO 10\n"
        "  ENDIF\n"
        "  IF (NG.GE.MXG) GO TO 20\n"
        "  NG = NG + 1\n"
        "  GSPCD(NG) = IFIA\n"
        "  TB(NG,1)=B0; TB(NG,2)=B1; TB(NG,3)=B2; TB(NG,4)=B3; TB(NG,5)=B4; TB(NG,6)=B5\n"
    )
    need(s, old_read, GOMPMORT)
    s = s.replace(old_read, new_read, 1)

    # -- GOMPLOAD: resolve onto FVS species index, copy 6th slot --
    old_res = ("        GB(ISPC,4)=TB(J,4); GB(ISPC,5)=TB(J,5)\n")
    need(s, old_res, GOMPMORT)
    s = s.replace(old_res,
                  "        GB(ISPC,4)=TB(J,4); GB(ISPC,5)=TB(J,5); GB(ISPC,6)=TB(J,6)\n", 1)

    # -- GOMPSURV: signature gains DBHV --
    need(s, "SUBROUTINE GOMPSURV(ISPC, CR, CCHV, FINTL, SURV)\n", GOMPMORT)
    s = s.replace("SUBROUTINE GOMPSURV(ISPC, CR, CCHV, FINTL, SURV)\n",
                  "SUBROUTINE GOMPSURV(ISPC, CR, CCHV, FINTL, SURV, DBHV)\n", 1)

    # -- GOMPSURV: declarations gain DBHV and B5 --
    need(s, "REAL CR, CCHV, FINTL, SURV\n", GOMPMORT)
    need(s, "REAL B0,B1,B2,B3,B4,CRC,CCHC,ETA,HZ,CTERM,SFLR\n", GOMPMORT)
    s = s.replace("REAL CR, CCHV, FINTL, SURV\n",
                  "REAL CR, CCHV, FINTL, SURV, DBHV\n", 1)
    s = s.replace("REAL B0,B1,B2,B3,B4,CRC,CCHC,ETA,HZ,CTERM,SFLR\n",
                  "REAL B0,B1,B2,B3,B4,B5,CRC,CCHC,ETA,HZ,CTERM,SFLR\n", 1)

    # -- GOMPSURV: load b5 --
    need(s, "B0=GB(ISPC,1); B1=GB(ISPC,2); B2=GB(ISPC,3); B3=GB(ISPC,4); B4=GB(ISPC,5)\n", GOMPMORT)
    s = s.replace(
        "B0=GB(ISPC,1); B1=GB(ISPC,2); B2=GB(ISPC,3); B3=GB(ISPC,4); B4=GB(ISPC,5)\n",
        "B0=GB(ISPC,1); B1=GB(ISPC,2); B2=GB(ISPC,3); B3=GB(ISPC,4); B4=GB(ISPC,5)\n"
        "B5=GB(ISPC,6)\n", 1)

    # -- GOMPSURV: add the ln(DBH) size term to eta (DBH in inches) --
    need(s, "ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM\n", GOMPMORT)
    s = s.replace(
        "ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM\n",
        "ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM + B5*LOG(MAX(DBHV,0.1))\n", 1)
    return s

# ---------------------------------------------------------------------------
# 4. vls/morts.f90 : pass tree diameter D (=DBH(I)) into GOMPSURV
#    (vls/morts.f90 is the mortality driver linked into FVSne)
# ---------------------------------------------------------------------------
VLSMORTS = "src-converted/vls/morts.f90"
def p_vlsmorts(s):
    if "CALL GOMPSURV(ISPC,CRG,CCHG,FINT,SURVG,D)" in s:
        return s
    need(s, "CALL GOMPSURV(ISPC,CRG,CCHG,FINT,SURVG)", VLSMORTS)
    s = s.replace("CALL GOMPSURV(ISPC,CRG,CCHG,FINT,SURVG)",
                  "CALL GOMPSURV(ISPC,CRG,CCHG,FINT,SURVG,D)", 1)
    return s

def main():
    for p in (GOMPMC, GOMPMORT, VLSMORTS):
        if not os.path.exists(p):
            sys.exit("ABORT: missing " + p)
    rw(GOMPMC, p_gompmc)
    rw(GOMPMORT, p_gompmort)
    rw(VLSMORTS, p_vlsmorts)
    print("DONE.")

if __name__ == "__main__":
    main()
