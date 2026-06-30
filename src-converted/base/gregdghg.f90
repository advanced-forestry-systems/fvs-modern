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
!     FVS_GREG_EMT      per-stand extreme min temperature (deg C).
!     FVS_GREG_TD       per-stand temperature difference MWMT-MCMT (deg C).
!     FVS_GREG_ELEV     stand elevation (feet); if unset, 0.
!==============================================================================
SUBROUTINE GREGLOADDG
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'GREGMC.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG)
REAL    TB(MXG,7)
CHARACTER(LEN=256) CVAL, CPATH
CHARACTER(LEN=512) LINE
INTEGER J, NG, IOS, U, IFIA, NN, ISPC
REAL B0,B1,B2,B3,B4,B5,B6
LOGICAL, SAVE :: LDONE = .FALSE.
!
IF (LDONE) RETURN
CALL GETENV('FVS_GREGDG', CVAL)
IF (CVAL.EQ.' ' .OR. CVAL(1:1).EQ.'0' .OR. CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N') THEN
  LGREGDG = .FALSE.; RETURN
ENDIF
LGREGDG = .TRUE.
LDONE = .TRUE.
!
GEMT = 0.0; GTD = 0.0; GELEV = 0.0
CALL GETENV('FVS_GREG_EMT', CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GEMT
CALL GETENV('FVS_GREG_TD',  CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GTD
CALL GETENV('FVS_GREG_ELEV',CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GELEV
!
CALL GETENV('FVS_GREGDG_COEF', CPATH)
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GREGDG: FVS_GREGDG set but FVS_GREGDG_COEF empty; NOT enabled.'
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
  READ(U,*,IOSTAT=IOS) IFIA, NN, B0, B1, B2, B3, B4, B5, B6
  IF (IOS.NE.0) GO TO 20
  IF (NG.GE.MXG) GO TO 20
  NG = NG + 1
  GSPCD(NG) = IFIA
  TB(NG,1)=B0; TB(NG,2)=B1; TB(NG,3)=B2; TB(NG,4)=B3
  TB(NG,5)=B4; TB(NG,6)=B5; TB(NG,7)=B6
  GO TO 10
20 CONTINUE
CLOSE(U)
NGREGDG = NG
!
DO ISPC=1,MAXSP
  GHAVE_DG(ISPC) = .FALSE.
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
        GHAVE_DG(ISPC) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
WRITE(JOSTND,*) 'GREGDG enabled: ', NG, ' species; EMT=', GEMT, ' TD=', GTD, ' ELEV=', GELEV
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
RETURN
END
