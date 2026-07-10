!==============================================================================
!  greghghg.f90 -- Greg Johnson four-arm height-growth + crown recession hooks.
!
!  GREGLOADHG loads per-species HG coefficients and per-stand climate, with arm
!  selection via the HGDRIVER keyword (IHGDRV in /GREGKW/). GREGHGV evaluates
!  the Chapman-Richards ORGANON annual height increment. GREGLOADCRW loads the
!  3-param crown recession coefficients. GREGCRWV evaluates the annual change in
!  height to live crown.
!
!  HG activation:
!    HGDRIVER keyword (IHGDRV): self-contained; selects CSV from confdir.
!      1 -> greg_hg_coefficients_bayes.csv       (Bayesian regional variants)
!      2 -> greg_hg_coefficients.csv             (Greg ORGANON CONUS -- default)
!      3 -> greg_hg_coefficients_conus_sp.csv    (CONUS species-dependent)
!      4 -> greg_hg_coefficients_conus_all.csv   (CONUS species-independent)
!    Env var fallback: FVS_GREGHG=1 + FVS_GREGHG_COEF=<path>
!
!  Crown recession activation:
!    CROWNDRIVER keyword (ICRWDRV=1): uses confdir/greg_crown_change_coefficients.csv
!    Env var fallback: FVS_GREGCRW=1 + FVS_GREGCRW_COEF=<path>
!
!  Shared climate env vars (read by GREGLOADHG; also read by GREGLOADDG for DG):
!    FVS_GREG_EMT      extreme min temperature (deg C)
!    FVS_GREG_TD       temperature diff MWMT-MCMT (deg C) -- shared via GTD in /GREGMR/
!    FVS_GREG_ELEV     stand elevation (feet)
!    FVS_GREG_CONFDIR  directory holding all config CSVs; default shown below
!
!  HG equation (Chapman-Richards ORGANON form, Greg Johnson 2026):
!    dHT_annual = MX * B1 * B2 * CR^B3
!                 * exp(-B1*HT - B4*CCFL - B8*sqrt(CCH) - B5*ELEV
!                       + B6*sqrt(TD) + B7*EMT)
!                 * (1 - exp(-B1*HT))^(B2-1)
!    GHG(ISPC,1)=MX, (2)=B1, (3)=B2, (4)=B3, (5)=B4,
!                (6)=B5, (7)=B6, (8)=B7, (9)=B8
!
!  Crown recession equation (Johnson, Marshall, Weiskittel 2026):
!    delta_htlc = crown_length * (1 - exp(b0 + b1*dHT_annual + b2*cch))
!    GCRW(ISPC,1)=b0, (2)=b1, (3)=b2
!==============================================================================

SUBROUTINE GREGLOADHG
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'GREGMC.f90'
INCLUDE 'GREGKW.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG)
REAL    TB(MXG,9)
CHARACTER(LEN=256) CVAL, CPATH, CDIR
CHARACTER(LEN=512) LINE
INTEGER J, NG, IOS, U, IFIA, NN, ISPC, JARM
REAL C0,C1,C2,C3,C4,C5,C6,C7,C8
LOGICAL LKWSEL
LOGICAL, SAVE :: LDONE = .FALSE.
!
IF (LDONE) RETURN
!
!  -- Enable check: HGDRIVER keyword (IHGDRV 1-4) or FVS_GREGHG env var --
LKWSEL = (IHGDRV.GE.1 .AND. IHGDRV.LE.4)
CALL GETENV('FVS_GREGHG', CVAL)
IF (.NOT.LKWSEL) THEN
  IF (CVAL.EQ.' ' .OR. CVAL(1:1).EQ.'0' .OR. &
      CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N') THEN
    LGREGHG = .FALSE.; RETURN
  ENDIF
ENDIF
LGREGHG = .TRUE.
LDONE   = .TRUE.
!
!  -- Climate scalars --
GEMT = 0.0; GELEV = 0.0
CALL GETENV('FVS_GREG_EMT',  CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GEMT
CALL GETENV('FVS_GREG_TD',   CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GTD
CALL GETENV('FVS_GREG_ELEV', CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GELEV
!
!  -- Resolve config directory --
CALL GETENV('FVS_GREG_CONFDIR', CDIR)
IF (CDIR.EQ.' ') CDIR = '/users/PUOM0008/crsfaaron/wt-dgdriver/config'
!
!  -- Resolve coefficient file path --
CPATH = ' '
IF (LKWSEL) THEN
  JARM = IHGDRV
  IF (JARM.EQ.1) CPATH = TRIM(CDIR)//'/greg_hg_coefficients_bayes.csv'
  IF (JARM.EQ.2) CPATH = TRIM(CDIR)//'/greg_hg_coefficients.csv'
  IF (JARM.EQ.3) CPATH = TRIM(CDIR)//'/greg_hg_coefficients_conus_sp.csv'
  IF (JARM.EQ.4) CPATH = TRIM(CDIR)//'/greg_hg_coefficients_conus_all.csv'
  WRITE(JOSTND,*) 'GREGHG: HGDRIVER arm ', JARM, ' -> ', TRIM(CPATH)
ELSE
  CALL GETENV('FVS_GREGHG_COEF', CPATH)
ENDIF
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GREGHG: no coefficient file (HGDRIVER unset and FVS_GREGHG_COEF empty); NOT enabled.'
  LGREGHG = .FALSE.; RETURN
ENDIF
!
!  -- Read HG coefficient CSV (SPCD, n, B0=MX, B1..B8) --
U = 69
OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
  WRITE(JOSTND,*) 'GREGHG: cannot open ', TRIM(CPATH), '; NOT enabled.'
  LGREGHG = .FALSE.; RETURN
ENDIF
READ(U,'(A)',IOSTAT=IOS) LINE      ! skip header
NG = 0
110 CONTINUE
  READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2, C3, C4, C5, C6, C7, C8
  IF (IOS.NE.0) GO TO 120
  IF (NG.GE.MXG) GO TO 120
  NG = NG + 1
  GSPCD(NG) = IFIA
  TB(NG,1)=C0; TB(NG,2)=C1; TB(NG,3)=C2; TB(NG,4)=C3; TB(NG,5)=C4
  TB(NG,6)=C5; TB(NG,7)=C6; TB(NG,8)=C7; TB(NG,9)=C8
  GO TO 110
120 CONTINUE
CLOSE(U)
NGREGHG = NG
!
!  -- Map FIAJSP -> GHG per FVS species slot --
DO ISPC=1,MAXSP
  GHAVE_HG(ISPC) = .FALSE.
  IFIA = -1
  IF (FIAJSP(ISPC).NE.' ') THEN
    READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
    IF (IOS.NE.0) IFIA = -1
  ENDIF
  IF (IFIA.GT.0) THEN
    DO J=1,NG
      IF (GSPCD(J).EQ.IFIA) THEN
        GHG(ISPC,1)=TB(J,1); GHG(ISPC,2)=TB(J,2); GHG(ISPC,3)=TB(J,3)
        GHG(ISPC,4)=TB(J,4); GHG(ISPC,5)=TB(J,5); GHG(ISPC,6)=TB(J,6)
        GHG(ISPC,7)=TB(J,7); GHG(ISPC,8)=TB(J,8); GHG(ISPC,9)=TB(J,9)
        GHAVE_HG(ISPC) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
!
!  -- Arm 4 override: species-independent -- apply single row to all species --
IF (IHGDRV.EQ.4 .AND. NG.GE.1) THEN
  DO ISPC=1,MAXSP
    GHG(ISPC,1)=TB(1,1); GHG(ISPC,2)=TB(1,2); GHG(ISPC,3)=TB(1,3)
    GHG(ISPC,4)=TB(1,4); GHG(ISPC,5)=TB(1,5); GHG(ISPC,6)=TB(1,6)
    GHG(ISPC,7)=TB(1,7); GHG(ISPC,8)=TB(1,8); GHG(ISPC,9)=TB(1,9)
    GHAVE_HG(ISPC) = .TRUE.
  ENDDO
ENDIF
!
!  -- Load optional MCW table for CCFL computation --
DO ISPC=1,MAXSP
  GHAVE_MCW(ISPC) = .FALSE.
  GMCW(ISPC,1)=0.0; GMCW(ISPC,2)=0.0; GMCW(ISPC,3)=0.0; GMCW(ISPC,4)=0.0
ENDDO
CALL GETENV('FVS_GREG_MCW_COEF', CPATH)
IF (CPATH.EQ.' ' .AND. LKWSEL) THEN
  CPATH = TRIM(CDIR)//'/greg_mcw_coefficients.csv'
ENDIF
IF (CPATH.NE.' ') THEN
  U = 70
  OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
  IF (IOS.EQ.0) THEN
    READ(U,'(A)',IOSTAT=IOS) LINE
    NG = 0
130 CONTINUE
      READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2
      IF (IOS.NE.0) GO TO 140
      IF (NG.GE.MXG) GO TO 140
      NG = NG + 1
      GSPCD(NG) = IFIA
      TB(NG,1)=REAL(NN); TB(NG,2)=C0; TB(NG,3)=C1; TB(NG,4)=C2
      GO TO 130
140 CONTINUE
    CLOSE(U)
    DO ISPC=1,MAXSP
      IFIA = -1
      IF (FIAJSP(ISPC).NE.' ') THEN
        READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
        IF (IOS.NE.0) IFIA = -1
      ENDIF
      IF (IFIA.GT.0) THEN
        DO J=1,NG
          IF (GSPCD(J).EQ.IFIA) THEN
            GMCW(ISPC,1)=TB(J,1); GMCW(ISPC,2)=TB(J,2)
            GMCW(ISPC,3)=TB(J,3); GMCW(ISPC,4)=TB(J,4)
            GHAVE_MCW(ISPC) = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    WRITE(JOSTND,*) 'GREGHG: loaded MCW for ', NG, ' species.'
  ELSE
    WRITE(JOSTND,*) 'GREGHG: MCW file not found: ', TRIM(CPATH), '; using FVS native CCFL.'
  ENDIF
ENDIF
!
WRITE(JOSTND,*) 'GREGHG enabled: arm=', IHGDRV, ' EMT=', GEMT, &
     ' TD=', GTD, ' ELEV=', GELEV
RETURN
END


SUBROUTINE GREGHGV(ISPC, HT, CR, CCFL, CCH, HGV)
!  Annual height increment (ft/yr), Greg Chapman-Richards ORGANON form.
!  Climate (GEMT, GTD, GELEV) read from /GREGHR/ common loaded by GREGLOADHG.
!  Uses double precision internally to match validated reference.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GREGMC.f90'
INTEGER ISPC
REAL HT, CR, CCFL, CCH, HGV
REAL(KIND=8) :: MX,B1,B2,B3,B4,B5,B6,B7,B8,CRP,CCHP,TDP,ARG,HG
MX=GHG(ISPC,1); B1=GHG(ISPC,2); B2=GHG(ISPC,3); B3=GHG(ISPC,4)
B4=GHG(ISPC,5); B5=GHG(ISPC,6); B6=GHG(ISPC,7); B7=GHG(ISPC,8); B8=GHG(ISPC,9)
CRP  = REAL(CR,8);   IF (CRP .LT.1.0D-4) CRP  = 1.0D-4
CCHP = REAL(CCH,8);  IF (CCHP.LT.0.0D0)  CCHP = 0.0D0
TDP  = REAL(GTD,8);  IF (TDP .LT.0.0D0)  TDP  = 0.0D0
ARG  = -B1*REAL(HT,8) - B4*REAL(CCFL,8) - B8*CCHP**0.5D0 &
       -B5*REAL(GELEV,8) + B6*SQRT(TDP) + B7*REAL(GEMT,8)
HG   = MX*B1*B2*CRP**B3*EXP(ARG)*(1.0D0-EXP(-B1*REAL(HT,8)))**(B2-1.0D0)
IF (HG.LT.0.0D0) HG = 0.0D0
HGV  = REAL(HG)
RETURN
END


SUBROUTINE GREGLOADCRW
!  Load per-species crown recession coefficients (b0, b1, b2).
!  Enabled by CROWNDRIVER keyword (ICRWDRV=1) or env var FVS_GREGCRW=1.
!  CSV: greg_crown_change_coefficients.csv with header SPCD,n,b0,b1,b2.
!  Model applied per cycle in the crown hook: GREGCRWV.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'GREGMC.f90'
INCLUDE 'GREGKW.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG)
REAL    TB(MXG,3)
CHARACTER(LEN=256) CVAL, CPATH, CDIR
CHARACTER(LEN=512) LINE
INTEGER J, NG, IOS, U, IFIA, NN, ISPC
REAL C0,C1,C2
LOGICAL LKWSEL
LOGICAL, SAVE :: LDONE = .FALSE.
!
IF (LDONE) RETURN
LKWSEL = (ICRWDRV.EQ.1)
CALL GETENV('FVS_GREGCRW', CVAL)
IF (.NOT.LKWSEL) THEN
  IF (CVAL.EQ.' ' .OR. CVAL(1:1).EQ.'0' .OR. &
      CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N') THEN
    LGREG_CRW = .FALSE.; RETURN
  ENDIF
ENDIF
LGREG_CRW = .TRUE.
LDONE      = .TRUE.
!
CALL GETENV('FVS_GREG_CONFDIR', CDIR)
IF (CDIR.EQ.' ') CDIR = '/users/PUOM0008/crsfaaron/wt-dgdriver/config'
!
CPATH = ' '
IF (LKWSEL) THEN
  CPATH = TRIM(CDIR)//'/greg_crown_change_coefficients.csv'
  WRITE(JOSTND,*) 'GREGCRW: CROWNDRIVER 1 -> ', TRIM(CPATH)
ELSE
  CALL GETENV('FVS_GREGCRW_COEF', CPATH)
ENDIF
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GREGCRW: no coefficient file; NOT enabled.'
  LGREG_CRW = .FALSE.; RETURN
ENDIF
!
U = 71
OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
  WRITE(JOSTND,*) 'GREGCRW: cannot open ', TRIM(CPATH), '; NOT enabled.'
  LGREG_CRW = .FALSE.; RETURN
ENDIF
READ(U,'(A)',IOSTAT=IOS) LINE    ! skip header
NG = 0
210 CONTINUE
  READ(U,*,IOSTAT=IOS) IFIA, NN, C0, C1, C2
  IF (IOS.NE.0) GO TO 220
  IF (NG.GE.MXG) GO TO 220
  NG = NG + 1
  GSPCD(NG) = IFIA
  TB(NG,1)=C0; TB(NG,2)=C1; TB(NG,3)=C2
  GO TO 210
220 CONTINUE
CLOSE(U)
NGREG_CRW = NG
!
DO ISPC=1,MAXSP
  GHAVE_CRW(ISPC) = .FALSE.
  IFIA = -1
  IF (FIAJSP(ISPC).NE.' ') THEN
    READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
    IF (IOS.NE.0) IFIA = -1
  ENDIF
  IF (IFIA.GT.0) THEN
    DO J=1,NG
      IF (GSPCD(J).EQ.IFIA) THEN
        GCRW(ISPC,1)=TB(J,1); GCRW(ISPC,2)=TB(J,2); GCRW(ISPC,3)=TB(J,3)
        GHAVE_CRW(ISPC) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
WRITE(JOSTND,*) 'GREGCRW enabled: ', NG, ' species.'
RETURN
END


SUBROUTINE GREGCRWV(ISPC, CL, DHTANN, CCH, DHTLC)
!  Annual change in height to live crown (ft/yr), positive = crown recession.
!  Equation: delta_htlc = CL * (1 - exp(b0 + b1*dHT_annual + b2*cch))
!  Returns DHTLC=0 if species not fitted or CL<=0.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GREGMC.f90'
INTEGER ISPC
REAL CL, DHTANN, CCH, DHTLC
REAL ARG
DHTLC = 0.0
IF (.NOT.GHAVE_CRW(ISPC)) RETURN
IF (CL.LE.0.0) RETURN
ARG  = GCRW(ISPC,1) + GCRW(ISPC,2)*DHTANN + GCRW(ISPC,3)*CCH
DHTLC = CL * (1.0 - EXP(ARG))
IF (DHTLC.LT.0.0) DHTLC = 0.0  ! crown can recede but not descend via this model
RETURN
END


BLOCK DATA GREGKWBD
!  Initialize keyword-selected driver codes to -1 (unset) before any keyword run.
IMPLICIT NONE
INCLUDE 'GREGKW.f90'
DATA IDGDRV /-1/, IMORTDRV /-1/, IHGDRV /-1/, ICRWDRV /-1/
END
