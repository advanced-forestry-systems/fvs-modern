!==============================================================================
!  GOMPMORT -- Greg Johnson gompit mortality, in the FVS growth loop.
!
!  Substitutes Greg Johnson's CONUS gompit survival for FVS native mortality,
!  per tree, per cycle, INSIDE TREGRO->MORTS so growth and density interact with
!  the substituted mortality each cycle. Gompit owns mortality for species that
!  carry a fitted row; unfit species keep FVS native background mortality.
!
!  Model (per species), annual hazard with a cycle-length exposure:
!      eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4
!      S   = 1 - exp(-exp(eta))             (annual survival, gompit)
!      S_period = S ** FINT                 (compound over the cycle)
!
!  cch is crown closure at the subject tree's tip. CHANGED 2026-07-08: this is
!  now a DIRECT port of Greg Johnson's authoritative biometrics.utilities
!  compute_cch (github.com/gregjohnsonbiometrics/biometrics.utilities src/cch.cpp).
!  cch is crown area as a FRACTION of an acre (AREACON = 0.25*PI/43560, NO x100)
!  in the horizontal plane tangent to each subject tree's tip, summed over all
!  TALLER trees. The prior implementation used a CCF-scale (x100) 40-layer
!  interpolation table plus an affine map (cch = 0.062 + 0.0036*cch_hat); that
!  extrapolated past the fitted 0.3-1.1 range on dense even-aged cohorts and
!  drove b3*cch^b4 to total mortality (the collapse). This port removes the x100,
!  the table, and the affine map so CCHT(I) is on Greg's fitted scale directly.
!
!  Activation (no recompile to toggle): environment variables read once in
!  GOMPLOAD: FVS_GOMPIT (1/on/true) and FVS_GOMPIT_COEF (coefficient csv path).
!==============================================================================

SUBROUTINE GOMPLOAD
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'PLOT.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'GOMPMC.f90'
!
INTEGER, PARAMETER :: MXG = 600
INTEGER GSPCD(MXG)
REAL    TB(MXG,5)
CHARACTER(LEN=256) CVAL, CPATH
CHARACTER(LEN=512) LINE
INTEGER I, J, NG, IOS, U, IFIA, NN, ISPC
REAL B0,B1,B2,B3,B4
LOGICAL, SAVE :: LDONE = .FALSE.
LOGICAL LENABLE
!
IF (LDONE) RETURN
LENABLE = .TRUE.
CALL GETENV('FVS_GOMPIT', CVAL)
IF (CVAL.EQ.' ') LENABLE = .FALSE.
IF (CVAL(1:1).EQ.'0' .OR. CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N' &
    .OR. CVAL(1:1).EQ.'f' .OR. CVAL(1:1).EQ.'F') LENABLE = .FALSE.
IF (.NOT.(LENABLE .OR. LGOMPKW)) RETURN
LDONE = .TRUE.
!
LGOMP = .FALSE.
NGOMP = 0
DO ISPC=1,MAXSP
  GHAVE(ISPC) = .FALSE.
  GGRP(ISPC)  = 16
  DO J=1,5
    GB(ISPC,J) = 0.0
  ENDDO
ENDDO
!
CALL GETENV('FVS_GOMPIT_COEF', CPATH)
IF (CPATH.EQ.' ') THEN
  WRITE(JOSTND,*) 'GOMPMORT: FVS_GOMPIT set but FVS_GOMPIT_COEF empty;', &
       ' gompit mortality NOT enabled.'
  RETURN
ENDIF
!
U = 67
OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
  WRITE(JOSTND,*) 'GOMPMORT: cannot open coeff file ', TRIM(CPATH), &
       '; gompit mortality NOT enabled.'
  RETURN
ENDIF
READ(U,'(A)',IOSTAT=IOS) LINE          ! header
NG = 0
10 CONTINUE
  READ(U,*,IOSTAT=IOS) IFIA, NN, B0, B1, B2, B3, B4
  IF (IOS.NE.0) GO TO 20
  IF (NG.GE.MXG) GO TO 20
  NG = NG + 1
  GSPCD(NG) = IFIA
  TB(NG,1)=B0; TB(NG,2)=B1; TB(NG,3)=B2; TB(NG,4)=B3; TB(NG,5)=B4
  GO TO 10
20 CONTINUE
CLOSE(U)
!
! resolve onto FVS species index via FIAJSP; assign ORGANON group by FIA code
DO ISPC=1,MAXSP
  IFIA = -1
  IF (FIAJSP(ISPC).NE.' ') THEN
    READ(FIAJSP(ISPC),*,IOSTAT=IOS) IFIA
    IF (IOS.NE.0) IFIA = -1
  ENDIF
  IF (IFIA.GT.0) THEN
    IF (IFIA.LT.300) THEN
      GGRP(ISPC) = 1
    ELSE
      GGRP(ISPC) = 16
    ENDIF
    DO J=1,NG
      IF (GSPCD(J).EQ.IFIA) THEN
        GB(ISPC,1)=TB(J,1); GB(ISPC,2)=TB(J,2); GB(ISPC,3)=TB(J,3)
        GB(ISPC,4)=TB(J,4); GB(ISPC,5)=TB(J,5)
        GHAVE(ISPC) = .TRUE.
        NGOMP = NGOMP + 1
        GO TO 30
      ENDIF
    ENDDO
30  CONTINUE
  ENDIF
ENDDO
!
LGOMP = .TRUE.
WRITE(JOSTND,*) 'GOMPMORT enabled: ', NG, ' fitted species read, ', &
     NGOMP, ' matched to this variant.'
RETURN
END


SUBROUTINE GOMPSURV(ISPC, CR, CCHV, FINTL, SURV)
!  Period survival for one tree of FVS species ISPC. Caller guarantees
!  GHAVE(ISPC). Returns SURV in (0,1]. CCHV is now Greg-scale cch (fraction
!  of an acre) from the compute_cch port in GOMPCCH.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
INTEGER ISPC
REAL CR, CCHV, FINTL, SURV
REAL B0,B1,B2,B3,B4,CRC,CCHC,ETA,HZ,CTERM
B0=GB(ISPC,1); B1=GB(ISPC,2); B2=GB(ISPC,3); B3=GB(ISPC,4); B4=GB(ISPC,5)
CRC = CR
IF (CRC.LT.1.0E-4) CRC = 1.0E-4
IF (CRC.GT.1.0)    CRC = 1.0
CCHC = CCHV
IF (CCHC.LT.0.0) CCHC = 0.0
IF (CCHC.GT.0.0) THEN
  CTERM = CCHC**B4
ELSE
  CTERM = 0.0
ENDIF
ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM
IF (ETA.GT.30.0)  ETA = 30.0
IF (ETA.LT.-30.0) ETA = -30.0
HZ = 1.0 - EXP(-EXP(ETA))   ! Greg gompit: annual SURVIVAL
SURV = MAX(0.0,MIN(1.0,HZ)) ** FINTL   ! compound annual survival over the cycle
RETURN
END


SUBROUTINE GOMPCCH
!  Recompute per-tree crown closure at tip (CCHT) for the current cycle.
!  DIRECT port of Greg Johnson's biometrics.utilities compute_cch (src/cch.cpp):
!  for each subject tree, cch = sum over all TALLER trees of the crown area at
!  the plane of the subject's tip, expressed as a FRACTION of an acre
!  (AREACON = 0.25*PI/43560, NO x100). No 40-layer interpolation table and no
!  affine map: CCHT(I) is Greg-scale and feeds GOMPSURV directly. Crown geometry
!  (MCW -> LCW, height-to-largest-crown-width HL, crown width above) is the same
!  ORGANON SWO port used previously; only the accumulation scale and method
!  changed. O(n^2) over the treelist, fine for per-plot n.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'ARRAYS.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'GOMPMC.f90'
!
REAL, PARAMETER :: PIC = 3.14159265358979
REAL, PARAMETER :: AREACON = 0.25*PIC/43560.0   ! crown area -> fraction of an acre
! ORGANON SWO (version 1) crown-width parameters, groups 1..18
REAL MCWB0(18),MCWB1(18),MCWB2(18),MCWPK(18)
REAL LCWB1(18),LCWB2(18),LCWB3(18)
REAL CWAB1(18),CWAB2(18),CWAB3(18),DACB(18)
! per-tree cached crown geometry
REAL    VLCW(MAXTRE), VHL(MAXTRE), VHT(MAXTRE), VH2D(MAXTRE)
INTEGER VG(MAXTRE)
LOGICAL VOK(MAXTRE)
INTEGER I, J, G
REAL DBHV, HTV, CRV, CL, HCB, MCWV, LCWV, HL, DD
REAL H, RP, CW, ALPHA, CCHI
!
DATA MCWB0/4.6366,6.1880,3.4835,4.6600546,3.2837,4.5652,4.0,4.5652, &
  3.4298629,2.9793895,4.4443,4.4443,4.0953,3.0785639,3.3625,8.0, &
  2.9793895,2.9793895/
DATA MCWB1/1.6078,1.0069,1.343,1.0701859,1.2031,1.4147,1.65,1.4147, &
  1.3532302,1.5512443,1.7040,1.7040,2.3849,1.9242211,2.0303,1.53, &
  1.5512443,1.5512443/
DATA MCWB2/-0.009625,0.0,-0.0082544,0.0,-0.0071858,0.0,0.0,0.0,0.0, &
  -0.01416129,0.0,0.0,-0.011630,0.0,-0.0073307,0.0,-0.01416129, &
  -0.01416129/
DATA MCWPK/88.52,999.99,81.35,999.99,83.71,999.99,999.99,999.99, &
  999.99,54.77,999.99,999.99,102.53,999.99,138.93,999.99,54.77,54.77/
DATA LCWB1/0.0,0.0,0.355532,0.0,-0.251389,0.0,-0.251389,0.0,0.118621, &
  0.0,0.0,0.0,0.0,0.364811,0.0,0.3227140,0.0,0.0/
DATA LCWB2/0.00371834,0.00308402,0.0,0.00339675,0.00692512,0.0, &
  0.00692512,0.0,0.00384872,0.0,0.0111972,0.0207676,0.0,0.0,0.0,0.0, &
  0.0,0.0/
DATA LCWB3/0.808121,0.0,0.0,0.532418,0.985922,0.0,0.985922,0.0,0.0, &
  1.161440,0.0,0.0,1.47018,0.0,1.27196,0.0,1.161440,1.161440/
DATA CWAB1/0.929973,0.999291,0.755583,0.755583,0.629785,0.629785, &
  0.629785,0.629785,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5/
DATA CWAB2/-0.135212,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
  0.0,0.0,0.0,0.0,0.0/
DATA CWAB3/-0.0157579,-0.0314603,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
  0.0,0.0,0.0,0.0,0.0,0.0,0.0/
DATA DACB/0.062,0.028454,0.05,0.05,0.20,0.209806,0.20,0.209806, &
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
!
! Phase 1: cache per-tree crown geometry (LCW, height-to-largest-crown-width HL)
DO I=1,ITRN
  CCHT(I) = 0.0
  VOK(I)  = .FALSE.
  DBHV = DBH(I)
  HTV  = HT(I)
  CRV  = REAL(ICR(I))/100.0
  IF (HTV.LE.0.0 .OR. DBHV.LE.0.0) CYCLE
  IF (CRV.LE.0.0 .OR. CRV.GT.1.0)  CYCLE
  G = GGRP(ISP(I))
  IF (G.LT.1 .OR. G.GT.18) G = 16
  CL  = CRV*HTV
  HCB = HTV - CL
  DD = DBHV
  IF (DD.GT.MCWPK(G)) DD = MCWPK(G)
  IF (HTV.LT.4.501) THEN
    MCWV = HTV/4.5*MCWB0(G)
  ELSE
    MCWV = MCWB0(G) + MCWB1(G)*DD + MCWB2(G)*DD*DD
  ENDIF
  LCWV = MCWV * CRV**(LCWB1(G) + LCWB2(G)*CL + LCWB3(G)*(DBHV/HTV))
  HL = HTV - (1.0-DACB(G))*CRV*HTV        ! height to largest crown width (hlcw)
  VLCW(I) = LCWV
  VHL(I)  = HL
  VHT(I)  = HTV
  VH2D(I) = HTV/DBHV
  VG(I)   = G
  VOK(I)  = .TRUE.
ENDDO
!
! Phase 2: for each subject tree, sum crown area of all TALLER trees at its tip
DO I=1,ITRN
  IF (.NOT.VOK(I)) THEN
    CCHT(I) = 0.0
    CYCLE
  ENDIF
  H    = VHT(I)          ! subject tip height
  CCHI = 0.0
  DO J=1,ITRN
    IF (.NOT.VOK(J)) CYCLE
    IF (VHT(J).LE.H) CYCLE            ! only trees strictly taller than the tip
    IF (H.LE.VHL(J)) THEN
      CW = VLCW(J)                    ! plane is at or below largest crown width
    ELSE
      RP = (VHT(J) - H)/(VHT(J) - VHL(J))   ! relative position above lcw
      IF (RP.LE.0.0) THEN
        CW = 0.0
      ELSE
        ALPHA = CWAB1(VG(J)) + CWAB2(VG(J))*SQRT(RP) + CWAB3(VG(J))*VH2D(J)
        IF (ALPHA.GT.0.0) THEN
          CW = VLCW(J) * RP**ALPHA
        ELSE
          CW = 1.0                    ! matches Greg's _cwa alpha<=0 guard
        ENDIF
      ENDIF
    ENDIF
    CCHI = CCHI + (CW*CW)*(AREACON*PROB(J))
  ENDDO
  CCHT(I) = CCHI
ENDDO
RETURN
END


SUBROUTINE GOMPON
!  GOMPMORT keyword hook (called from vbase/initre.f90 when the GOMPMORT keyword
!  is read). Flags keyword activation; GOMPLOAD then loads coefficients.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
LGOMPKW = .TRUE.
RETURN
END


BLOCK DATA GOMPBD
!  Initialise the gompit COMMON flags so LGOMP / LGOMPKW / GHAVE are defined
!  before any activation.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
DATA LGOMP /.FALSE./, LGOMPKW /.FALSE./, GHAVE /MAXSP*.FALSE./
DATA NGOMP /0/, GGRP /MAXSP*16/
END
