!==============================================================================
!  gregdghg.f90 -- Greg Johnson native diameter-growth substitution.
!
!  GREGLOADDG loads the per-species Greg DG coefficients (B0..B6) and the per-stand
!  climate scalars, mirroring GOMPLOAD. GREGDGV evaluates the deployed annual
!  diameter increment (validated to 8e-8 vs the projector dg_annual). The DGF hook
!  (separate edit) calls GREGDGV over the cycle for covered species and rebuilds
!  the engine WK2 = log(inside-bark DDS) via BRATIO.
!
!  Activation (no recompile to toggle): env vars read once in GREGLOADDG:
!     FVS_GREGDG        1/on/true to enable Greg DG substitution.
!     FVS_GREGDG_COEF   path to greg_dg_coefficients.csv (header + SPCD,n,B0..B6).
!     FVS_GREG_CONFDIR  dir holding the config CSVs; used to resolve DGDRIVER-selected
!                       files. Default '/users/PUOM0008/crsfaaron/wt-dgdriver/config'.
!     FVS_GREG_EMT      per-stand extreme min temperature (deg C).
!     FVS_GREG_TD       per-stand temperature difference MWMT-MCMT (deg C).
!     FVS_GREG_ELEV     stand elevation (feet); if unset, 0.
!
!  DGDRIVER keyword (common /GREGKW/ IDGDRV, GREGKW.f90) selects the coefficient
!  file directly, so a keyword-only run needs no FVS_GREGDG / FVS_GREGDG_COEF:
!     DGDRIVER 1 -> <confdir>/greg_dg_coefficients.csv       (deployed)
!     DGDRIVER 2 -> <confdir>/greg_dg_coefficients_refit.csv (refit on our data)
!  A mapped IDGDRV both turns the hook on (LGREGDG) and picks the file; it takes
!  precedence over FVS_GREGDG_COEF. Unmapped/unset IDGDRV falls back to the env
!  var behaviour. NOTE: the 8-column driver-family files (cspi/bgi/esi/elev/emt)
!  are NOT loadable by this 9-column reader, so codes other than 1/2 are not yet
!  wired here (downstream / out of scope).
!==============================================================================
SUBROUTINE GREGLOADDG
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'GREGMC.f90'
INCLUDE 'GREGKW.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG)
REAL    TB(MXG,7), TBMAX(MXG)
CHARACTER(LEN=256) CVAL, CPATH, CDIR
LOGICAL LKWSEL
CHARACTER(LEN=512) LINE
INTEGER J, NG, IOS, U, IFIA, NN, ISPC
REAL B0,B1,B2,B3,B4,B5,B6,BMX
REAL, PARAMETER :: DGMAX_NONE = 1.0E6   ! sentinel: no DBHMAX -> decel==1
LOGICAL, SAVE :: LDONE = .FALSE.
!
IF (LDONE) RETURN
!  LKWSEL: did the DGDRIVER keyword set a code we map to a loadable coef file?
!  Codes 1 (deployed) and 2 (refit) are the two engine-loadable 9-col files.
LKWSEL = (IDGDRV.EQ.1 .OR. IDGDRV.EQ.2)
CALL GETENV('FVS_GREGDG', CVAL)
!  Enable if either the env toggle is on OR a mapped DGDRIVER code was given.
!  This makes DGDRIVER self-contained: a keyword-only run needs no env var.
IF (.NOT.LKWSEL) THEN
  IF (CVAL.EQ.' ' .OR. CVAL(1:1).EQ.'0' .OR. CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N') THEN
    LGREGDG = .FALSE.; RETURN
  ENDIF
ENDIF
LGREGDG = .TRUE.
LDONE = .TRUE.
!
GEMT = 0.0; GTD = 0.0; GELEV = 0.0
CALL GETENV('FVS_GREG_EMT', CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GEMT
CALL GETENV('FVS_GREG_TD',  CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GTD
CALL GETENV('FVS_GREG_ELEV',CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GELEV
!
!  ---- Coefficient-file resolution ------------------------------------------
!  Precedence: a mapped DGDRIVER code selects the file (keyword self-contained);
!  otherwise fall back to FVS_GREGDG_COEF (the original env-var behaviour).
CPATH = ' '
IF (LKWSEL) THEN
  !  Resolve config dir robustly: FVS_GREG_CONFDIR if set, else a known abs dir.
  CALL GETENV('FVS_GREG_CONFDIR', CDIR)
  IF (CDIR.EQ.' ') CDIR = '/users/PUOM0008/crsfaaron/wt-dgdriver/config'
  !  DGDRIVER code -> coefficient filename mapping.
  IF (IDGDRV.EQ.1) THEN
    CPATH = TRIM(CDIR)//'/greg_dg_coefficients.csv'        ! deployed
  ELSE IF (IDGDRV.EQ.2) THEN
    CPATH = TRIM(CDIR)//'/greg_dg_coefficients_refit.csv'  ! refit on our data
  ENDIF
  WRITE(JOSTND,*) 'GREGDG: DGDRIVER code ', IDGDRV, ' selected ', TRIM(CPATH)
ELSE
  CALL GETENV('FVS_GREGDG_COEF', CPATH)
ENDIF
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GREGDG: no coefficient file (DGDRIVER unset/unmapped and ', &
       'FVS_GREGDG_COEF empty); NOT enabled.'
  LGREGDG = .FALSE.; RETURN
ENDIF
U = 68
OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
  WRITE(JOSTND,*) 'GREGDG: cannot open ', TRIM(CPATH), '; NOT enabled.'
  LGREGDG = .FALSE.; RETURN
ENDIF
READ(U,'(A)',IOSTAT=IOS) LINE                ! header
NG = 0
10 CONTINUE
  READ(U,'(A)',IOSTAT=IOS) LINE
  IF (IOS.NE.0) GO TO 20
  IF (LINE.EQ.' ') GO TO 10
  IF (NG.GE.MXG) GO TO 20
  !  Try 10-field record (SPCD,n,B0..B6,DBHMAX); if that fails, fall back to the
  !  legacy 9-field record and leave DBHMAX at the sentinel (decel==1).
  BMX = DGMAX_NONE
  READ(LINE,*,IOSTAT=IOS) IFIA, NN, B0, B1, B2, B3, B4, B5, B6, BMX
  IF (IOS.NE.0) THEN
    BMX = DGMAX_NONE
    READ(LINE,*,IOSTAT=IOS) IFIA, NN, B0, B1, B2, B3, B4, B5, B6
    IF (IOS.NE.0) GO TO 10
  ENDIF
  IF (BMX.LE.0.0) BMX = DGMAX_NONE
  NG = NG + 1
  GSPCD(NG) = IFIA
  TB(NG,1)=B0; TB(NG,2)=B1; TB(NG,3)=B2; TB(NG,4)=B3
  TB(NG,5)=B4; TB(NG,6)=B5; TB(NG,7)=B6
  TBMAX(NG)=BMX
  GO TO 10
20 CONTINUE
CLOSE(U)
NGREGDG = NG
!
DO ISPC=1,MAXSP
  GHAVE_DG(ISPC) = .FALSE.
  GDGMAX(ISPC) = DGMAX_NONE
  IFIA = -1
  IF (FIAJSP(ISPC).NE.' ') THEN
    READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
    IF (IOS.NE.0) IFIA = -1
  ENDIF
  IF (IFIA.GT.0) THEN
    DO J=1,NG
      IF (GSPCD(J).EQ.IFIA) THEN
        GDG(ISPC,1)=TB(J,1); GDG(ISPC,2)=TB(J,2); GDG(ISPC,3)=TB(J,3)
        GDG(ISPC,4)=TB(J,4); GDG(ISPC,5)=TB(J,5); GDG(ISPC,6)=TB(J,6); GDG(ISPC,7)=TB(J,7)
        GDGMAX(ISPC) = TBMAX(J)
        GHAVE_DG(ISPC) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
WRITE(JOSTND,*) 'GREGDG enabled: ', NG, ' species; EMT=', GEMT, ' TD=', GTD, ' ELEV=', GELEV
IF (IDGDRV.GE.0) WRITE(JOSTND,*) 'GREGDG keyword DGDRIVER code=', IDGDRV
RETURN
END

SUBROUTINE GREGDGV(ISPC, DBH, CR, HT, BAL, ELEV, EMT, G)
!  Annual diameter increment (inches/yr), Greg deployed DG. Validated 8e-8 vs the
!  projector dg_annual.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GREGMC.f90'
INTEGER ISPC
REAL DBH, CR, HT, BAL, ELEV, EMT, G, Z, CRC, HTC, BALC, ARGNUM, ARGDEN
REAL DECEL, DBHMX, XR
CRC = CR;  IF (CRC.LT.1.0E-4) CRC = 1.0E-4
HTC = HT;  IF (HTC.LT.0.0) HTC = 0.0
BALC = BAL; IF (BALC.LT.0.0) BALC = 0.0
ARGNUM = (DBH+1.0)**2
ARGDEN = (CRC*HTC+1.0)**GDG(ISPC,4)
Z = GDG(ISPC,1) + GDG(ISPC,2)*LOG(ARGNUM/ARGDEN) &
    + GDG(ISPC,3)*BALC**GDG(ISPC,5)/LOG(DBH+2.7) &
    + GDG(ISPC,6)*ELEV + GDG(ISPC,7)*EMT
IF (Z.GT.5.0)   Z = 5.0
IF (Z.LT.-30.0) Z = -30.0
G = EXP(Z); IF (G.LT.0.0) G = 0.0
!  ---- Size-based deceleration (COR-style plateau) --------------------------
!  Native NE-TWIGS DG plateaus via a size calibration this hook otherwise omits,
!  so an unbounded compounding loop runs QMD away over multi-century horizons.
!  Multiply the annual increment by a logistic that is ~1 until DBH nears the
!  per-species maximum diameter GDGMAX (from the DBHMAX coef column) and ramps
!  to 0 as DBH -> GDGMAX. Onset ~85% of max, half-width 4% of max, so trees far
!  below max (all short remeasurement intervals) are essentially untouched.
!  Sentinel GDGMAX (no DBHMAX column) leaves DECEL=1 -> exact old behaviour.
DBHMX = GDGMAX(ISPC)
IF (DBHMX .LT. 1.0E5) THEN
  XR = (DBH - 0.85*DBHMX) / (0.04*DBHMX)
  IF (XR .GT. 30.0) THEN
    DECEL = 0.0
  ELSE IF (XR .LT. -30.0) THEN
    DECEL = 1.0
  ELSE
    DECEL = 1.0 / (1.0 + EXP(XR))
  ENDIF
  IF (DECEL .LT. 0.0) DECEL = 0.0
  IF (DECEL .GT. 1.0) DECEL = 1.0
  G = G * DECEL
ENDIF
RETURN
END

BLOCK DATA GREGKWBD
!  Initialise keyword-selected driver codes to -1 (unset) before any keyword runs.
IMPLICIT NONE
INCLUDE 'GREGKW.f90'
DATA IDGDRV /-1/, IMORTDRV /-1/
END
