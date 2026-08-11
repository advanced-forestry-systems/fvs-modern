!==============================================================================
!  greghtdbh.f90 -- Marshall CONUS native height-diameter (HT-DBH) substitution.
!
!  GREGLOADHD loads the per-species CONUS HT-DBH coefficients (model_id, B1..B3),
!  mirroring GREGLOADDG. Function GREGHTDBH evaluates the predicted total height
!  (feet) from DBH (inches) for a covered species, dispatching one of the six
!  Marshall model forms. The HTDBH hook (edit in ne/htdbh.f90) calls GREGHTDBH in
!  MODE 0 (height-from-diameter) for covered species and substitutes its height
!  for the native Wykoff / Curtis-Arney prediction; uncovered species and other
!  modes are unchanged.
!
!  Activation (no recompile to toggle): env vars read once in GREGLOADHD:
!     FVS_GREGHTDBH        1/on/true to enable the CONUS HT-DBH substitution.
!     FVS_GREGHTDBH_COEF   path to conus_htdbh_coefficients.csv
!                          (header + SPCD,model_id,B1,B2,B3,...). B3 may be blank
!                          for models 1 and 4. Only the first five columns are read.
!  Default OFF and silent when FVS_GREGHTDBH is unset.
!
!  The six forms (bh = 4.5 ft; DBH in inches; HT returned in feet), transcribed
!  verbatim from the model.formula strings in the Marshall CONUS_HDfit_Model_*.CSV:
!     model 1: HT = bh + exp(B1 + B2*DBH**(-1.0))
!     model 2: HT = bh + exp(B1 + B2*DBH**B3)
!     model 3: HT = bh + B1*(1.0-exp(B2*DBH))**B3          (Chapman-Richards)
!     model 4: HT = bh + exp(B1 + B2/(DBH+1.0))
!     model 5: HT = bh + exp(B1 + B2/(DBH+B3))
!     model 6: HT = bh + B1*(1.0-exp(B2*DBH**B3))          (Weibull)
!  Species not covered => GREGHTDBH returns a sentinel (<=0) so the caller keeps
!  the native value.
!==============================================================================
SUBROUTINE GREGLOADHD
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'GREGHDMC.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG), GMOD(MXG)
REAL    TB(MXG,3)
CHARACTER(LEN=256) CVAL, CPATH
CHARACTER(LEN=1024) LINE
INTEGER J, NG, IOS, U, IFIA, IMOD, ISPC
REAL B1, B2, B3
LOGICAL, SAVE :: LDONE = .FALSE.
!
IF (LDONE) RETURN
CALL GETENV('FVS_GREGHTDBH', CVAL)
IF (CVAL.EQ.' ' .OR. CVAL(1:1).EQ.'0' .OR. CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N') THEN
  LGREGHTDBH = .FALSE.
  RETURN
ENDIF
LGREGHTDBH = .TRUE.
LDONE = .TRUE.
!
CALL GETENV('FVS_GREGHTDBH_COEF', CPATH)
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GREGHTDBH: FVS_GREGHTDBH_COEF empty; NOT enabled.'
  LGREGHTDBH = .FALSE.; RETURN
ENDIF
U = 69
OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
  WRITE(JOSTND,*) 'GREGHTDBH: cannot open ', TRIM(CPATH), '; NOT enabled.'
  LGREGHTDBH = .FALSE.; RETURN
ENDIF
READ(U,'(A)',IOSTAT=IOS) LINE                 ! header
NG = 0
10 CONTINUE
  READ(U,'(A)',IOSTAT=IOS) LINE
  IF (IOS.NE.0) GO TO 20
  IF (LINE.EQ.' ') GO TO 10
  IF (NG.GE.MXG) GO TO 20
  !  Parse only SPCD,model_id,B1,B2,B3 (first five columns). B3 is blank for
  !  models 1 and 4; the blank field parses to 0.0 (unused by those forms).
  CALL GREGHD_PARSE5(LINE, IFIA, IMOD, B1, B2, B3, IOS)
  IF (IOS.NE.0) GO TO 10
  IF (IMOD.LT.1 .OR. IMOD.GT.6) GO TO 10
  NG = NG + 1
  GSPCD(NG) = IFIA
  GMOD(NG)  = IMOD
  TB(NG,1)=B1; TB(NG,2)=B2; TB(NG,3)=B3
  GO TO 10
20 CONTINUE
CLOSE(U)
!
DO ISPC=1,MAXSP
  GHAVE_HD(ISPC) = .FALSE.
  GHDMOD(ISPC) = 0
  GHD(ISPC,1)=0.0; GHD(ISPC,2)=0.0; GHD(ISPC,3)=0.0
  IFIA = -1
  IF (FIAJSP(ISPC).NE.' ') THEN
    READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
    IF (IOS.NE.0) IFIA = -1
  ENDIF
  IF (IFIA.GT.0) THEN
    DO J=1,NG
      IF (GSPCD(J).EQ.IFIA) THEN
        GHDMOD(ISPC) = GMOD(J)
        GHD(ISPC,1)=TB(J,1); GHD(ISPC,2)=TB(J,2); GHD(ISPC,3)=TB(J,3)
        GHAVE_HD(ISPC) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
NGREGHD = NG
WRITE(JOSTND,*) 'GREGHTDBH enabled: ', NG, ' species loaded from ', TRIM(CPATH)
RETURN
END

SUBROUTINE GREGHD_PARSE5(LINE, IFIA, IMOD, B1, B2, B3, IOS)
!  Parse the first five comma-separated fields (SPCD,model_id,B1,B2,B3) of a
!  coefficient record, tolerating a blank B3 field (models 1 and 4) which becomes
!  0.0. Avoids list-directed READ swallowing across an empty field.
IMPLICIT NONE
CHARACTER(LEN=*) LINE
INTEGER IFIA, IMOD, IOS
REAL B1, B2, B3
INTEGER P0, P1, NF, LN
CHARACTER(LEN=256) FLD
REAL RVAL
INTEGER IVAL
IFIA=0; IMOD=0; B1=0.0; B2=0.0; B3=0.0; IOS=0
LN = LEN_TRIM(LINE)
P0 = 1
NF = 0
DO
  P1 = INDEX(LINE(P0:LN), ',')
  IF (P1.EQ.0) THEN
    FLD = LINE(P0:LN)
  ELSE
    FLD = LINE(P0:P0+P1-2)
  ENDIF
  NF = NF + 1
  IF (FLD.EQ.' ') THEN
    RVAL = 0.0; IVAL = 0
  ELSE
    IF (NF.LE.2) THEN
      READ(FLD,*,IOSTAT=IOS) IVAL
      IF (IOS.NE.0) RETURN
    ELSE
      READ(FLD,*,IOSTAT=IOS) RVAL
      IF (IOS.NE.0) RETURN
    ENDIF
  ENDIF
  SELECT CASE (NF)
    CASE (1); IFIA = IVAL
    CASE (2); IMOD = IVAL
    CASE (3); B1 = RVAL
    CASE (4); B2 = RVAL
    CASE (5); B3 = RVAL
  END SELECT
  IF (NF.GE.5) EXIT
  IF (P1.EQ.0) EXIT               ! ran out of fields before 5
  P0 = P0 + P1
  IF (P0.GT.LN) EXIT
ENDDO
IF (NF.LT.4) IOS = 1              ! need at least SPCD,model,B1,B2
RETURN
END

REAL FUNCTION GREGHTDBH(ISPC, DBH)
!  CONUS HT-DBH predicted total height (feet) from DBH (inches) for FVS species
!  ISPC. Dispatches one of the six Marshall model forms by GHDMOD(ISPC). Returns
!  a sentinel (<=0.0) for uncovered species so the caller keeps the native value.
!  Compiled real(4) here (engine default REAL); the standalone verification
!  harness re-runs the same forms in real(8) and clears the <1e-6 gate.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GREGHDMC.f90'
INTEGER ISPC
REAL DBH, BH, B1, B2, B3, DD, HT
BH = 4.5
GREGHTDBH = -1.0
IF (ISPC.LT.1 .OR. ISPC.GT.MAXSP) RETURN
IF (.NOT.GHAVE_HD(ISPC)) RETURN
DD = DBH
IF (DD.LE.0.0) RETURN
B1 = GHD(ISPC,1); B2 = GHD(ISPC,2); B3 = GHD(ISPC,3)
SELECT CASE (GHDMOD(ISPC))
  CASE (1)
    HT = BH + EXP(B1 + B2*DD**(-1.0))
  CASE (2)
    HT = BH + EXP(B1 + B2*DD**B3)
  CASE (3)
    HT = BH + B1*(1.0-EXP(B2*DD))**B3
  CASE (4)
    HT = BH + EXP(B1 + B2/(DD+1.0))
  CASE (5)
    HT = BH + EXP(B1 + B2/(DD+B3))
  CASE (6)
    HT = BH + B1*(1.0-EXP(B2*DD**B3))
  CASE DEFAULT
    RETURN
END SELECT
GREGHTDBH = HT
RETURN
END
