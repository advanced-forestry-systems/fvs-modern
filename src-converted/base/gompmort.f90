!==============================================================================
!  GOMPMORT -- Greg Johnson gompit mortality, in the FVS growth loop.
!
!  Substitutes Greg Johnson's CONUS gompit survival for FVS native mortality,
!  per tree, per cycle, INSIDE TREGRO->MORTS so growth and density interact with
!  the substituted mortality each cycle (unlike a post-hoc TPA overlay, which
!  cannot regulate density). Gompit owns mortality for species that carry a
!  fitted row; unfit species keep FVS native background mortality.
!
!  Model (per species), annual hazard with a cycle-length exposure:
!      eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4
!      H   = exp(eta)                       (annual hazard)
!      S   = exp(-H * FINT)                 (period survival, FINT = cycle yrs)
!      trees killed = PROB * (1 - S)
!
!  cch is crown closure at the subject tree's tip, recomputed each cycle from
!  the live treelist by an ORGANON crown-closure port (GOMPCCH), then mapped
!  onto the gompit cch scale by the validated affine fit
!      cch = CCH_A + CCH_B * cch_hat        (35d_validate_cch.R, Spearman 0.93).
!
!  AARON'S VALIDATED ASSUMPTIONS (science owned by Aaron Weiskittel; this is the
!  Fortran wiring of his model):
!    * gompit coefficients (greg_mortality_coefficients.csv).
!    * ORGANON SWO crown geometry as the cch proxy (CAL_CCH.for lineage).
!    * coarse FIA-species -> ORGANON group map: softwood (FIA<300) -> 1 (DF),
!      hardwood -> 16 (RA). Refine GGRP for a tighter fit.
!    * affine map CCH_A=0.062, CCH_B=0.0036.
!
!  Activation (no recompile to toggle): environment variables read once in
!  GOMPLOAD (called from MORCON):
!      FVS_GOMPIT       set to 1/on/true to enable.
!      FVS_GOMPIT_COEF  path to greg_mortality_coefficients.csv
!                       (columns: SPCD,n,b0,b1,b2,b3,b4,...).
!  A GOMPMORT keyword can be layered on later; the env switch keeps the first
!  integration testable and gives a clean default-vs-gompit A/B.
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
! resolve once per run (COMMON state persists across the per-stand MORCON calls).
! Activation comes from the FVS_GOMPIT env var OR the GOMPMORT keyword (LGOMPKW,
! set by GOMPON in initre, which runs before MORCON). Do not lock LDONE until
! actually activated, so a keyword seen on a later stand can still turn it on.
IF (LDONE) RETURN
LENABLE = .TRUE.
CALL GETENV('FVS_GOMPIT', CVAL)
IF (CVAL.EQ.' ') LENABLE = .FALSE.
IF (CVAL(1:1).EQ.'0' .OR. CVAL(1:1).EQ.'n' .OR. CVAL(1:1).EQ.'N' &
    .OR. CVAL(1:1).EQ.'f' .OR. CVAL(1:1).EQ.'F') LENABLE = .FALSE.
IF (.NOT.(LENABLE .OR. LGOMPKW)) RETURN
LDONE = .TRUE.
!
! -- Small-DBH mortality guard (young-cohort collapse fix). Read once here so
! -- the threshold/floor are tunable without a recompile. Defaults keep the
! -- guard ON but leave every tree with DBH>=GDBHMIN byte-identical to prior
! -- gompit. FVS_GOMP_DBHMIN<=0 disables the guard entirely.
GDBHMIN = 1.0
GSFLOOR = 0.95
CALL GETENV('FVS_GOMP_DBHMIN', CVAL)
IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GDBHMIN
CALL GETENV('FVS_GOMP_SFLOOR', CVAL)
IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GSFLOOR
IF (GSFLOOR.LT.0.0) GSFLOOR = 0.0
IF (GSFLOOR.GT.1.0) GSFLOOR = 1.0
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


SUBROUTINE GOMPSURV(ISPC, DBHV, CR, CCHV, FINTL, SURV)
!  Period survival for one tree of FVS species ISPC. Caller guarantees
!  GHAVE(ISPC). Returns SURV in (0,1].
!
!  SMALL-DBH GUARD (young-cohort collapse fix): the raw gompit annual survival
!  1-exp(-exp(eta)) can drop toward 0 for some conifers at the crown ratio and
!  crown closure of a young dense seedling cohort, killing ~100% of the stand
!  in the first cycle (falsified by the Bakuzis assessment on NE/SN conifers).
!  For trees with 0 <= DBHV < GDBHMIN the per-cycle survival is FLOORED at
!  GSFLOOR**FINTL (annual survival GSFLOOR), so a seedling cohort self-thins at
!  a bounded plausible rate instead of collapsing. The floor is a MAX, so where
!  gompit already predicts survival above the floor it is left untouched. Trees
!  with DBHV >= GDBHMIN, and any caller passing DBHV < 0 (unknown DBH), take the
!  identical eta/hazard path as before: established-stand behavior is
!  BYTE-IDENTICAL. GDBHMIN and GSFLOOR come from FVS_GOMP_DBHMIN and
!  FVS_GOMP_SFLOOR (defaults 1.0 in and 0.95 per yr); GDBHMIN<=0 disables.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
INTEGER ISPC
REAL DBHV, CR, CCHV, FINTL, SURV
REAL B0,B1,B2,B3,B4,CRC,CCHC,ETA,HZ,CTERM,SFLR
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
HZ = 1.0 - EXP(-EXP(ETA))   ! Greg gompit: annual SURVIVAL (high eta -> high survival)
SURV = MAX(0.0,MIN(1.0,HZ)) ** FINTL   ! compound annual survival over the cycle
!
! Small-DBH guard: floor per-cycle survival for sub-threshold trees only.
! GDBHMIN<=0 disables. DBHV<0 (caller signals no DBH) bypasses the guard so
! behavior is identical to the pre-guard code path.
IF (GDBHMIN.GT.0.0 .AND. DBHV.GE.0.0 .AND. DBHV.LT.GDBHMIN) THEN
  SFLR = GSFLOOR ** FINTL
  IF (SURV.LT.SFLR) SURV = SFLR
ENDIF
RETURN
END


SUBROUTINE GOMPCCH
!  Recompute per-tree crown closure at tip (CCHT) for the current cycle, an
!  ORGANON-crown port (faithful to cch_organon.py so the validated affine map
!  applies). Fills CCHT(I) = CCH_A + CCH_B * cch_hat for every live tree.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'ARRAYS.f90'
INCLUDE 'CONTRL.f90'
INCLUDE 'GOMPMC.f90'
!
REAL, PARAMETER :: CCH_A = 0.062, CCH_B = 0.0036
! ORGANON SWO (version 1) crown-width parameters, groups 1..18
REAL MCWB0(18),MCWB1(18),MCWB2(18),MCWPK(18)
REAL LCWB1(18),LCWB2(18),LCWB3(18)
REAL CWAB1(18),CWAB2(18),CWAB3(18),DACB(18)
REAL CCH(0:40)
INTEGER I, G, II
REAL DBHV, HTV, CRV, EXPN, MAXHT, CL, HCB, MCWV, LCWV, HL
REAL XL, BND, CW, TOP, XI, XXI, CCHHAT, DD, RP, HtoD
INTEGER IDX
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
DO I=1,ITRN
  CCHT(I) = 0.0
ENDDO
!
! tallest valid tree
MAXHT = 0.0
DO I=1,ITRN
  IF (HT(I).GT.MAXHT .AND. DBH(I).GT.0.0 .AND. ICR(I).GT.0) MAXHT=HT(I)
ENDDO
IF (MAXHT.LE.0.0) RETURN
!
DO II=0,40
  CCH(II) = 0.0
ENDDO
CCH(40) = MAXHT
!
DO I=1,ITRN
  DBHV = DBH(I)
  HTV  = HT(I)
  CRV  = REAL(ICR(I))/100.0
  EXPN = PROB(I)
  IF (HTV.LE.0.0 .OR. DBHV.LE.0.0) CYCLE
  IF (CRV.LE.0.0 .OR. CRV.GT.1.0)  CYCLE
  G = GGRP(ISP(I))
  IF (G.LT.1 .OR. G.GT.18) G = 16
  CL  = CRV*HTV
  HCB = HTV - CL
  ! MCW
  DD = DBHV
  IF (DD.GT.MCWPK(G)) DD = MCWPK(G)
  IF (HTV.LT.4.501) THEN
    MCWV = HTV/4.5*MCWB0(G)
  ELSE
    MCWV = MCWB0(G) + MCWB1(G)*DD + MCWB2(G)*DD*DD
  ENDIF
  ! LCW
  LCWV = MCWV * CRV**(LCWB1(G) + LCWB2(G)*CL + LCWB3(G)*(DBHV/HTV))
  ! height to largest crown width
  HL = HTV - (1.0-DACB(G))*CRV*HTV
  HtoD = HTV/DBHV
  DO II=40,1,-1
    XL = REAL(II-1)*(CCH(40)/40.0)
    BND = HCB
    IF (HL.GT.HCB) BND = HL
    IF (XL.LE.BND) THEN
      IF (HCB.LE.HL) THEN
        CW = LCWV
      ELSE
        ! cw_above at max(XL,HCB)
        RP = (HTV - MAX(XL,HCB))/(HTV - HL)
        IF (RP.LE.0.0) THEN
          CW = 0.0
        ELSE
          CW = LCWV * RP**(CWAB1(G)+CWAB2(G)*SQRT(RP)+CWAB3(G)*HtoD)
        ENDIF
      ENDIF
    ELSEIF (XL.LT.HTV) THEN
      RP = (HTV - XL)/(HTV - HL)
      IF (RP.LE.0.0) THEN
        CW = 0.0
      ELSE
        CW = LCWV * RP**(CWAB1(G)+CWAB2(G)*SQRT(RP)+CWAB3(G)*HtoD)
      ENDIF
    ELSE
      CW = 0.0
    ENDIF
    CCH(II) = CCH(II) + (CW*CW)*(0.001803*EXPN)
  ENDDO
ENDDO
!
! interpolate each tree's tip cch_hat, apply affine map
TOP = CCH(40)
DO I=1,ITRN
  HTV = HT(I)
  IF (DBH(I).LE.0.0 .OR. ICR(I).LE.0 .OR. HTV.LE.0.0) THEN
    CCHT(I) = CCH_A
    CYCLE
  ENDIF
  IF (HTV.GE.TOP .OR. TOP.LE.0.0) THEN
    CCHHAT = 0.0
  ELSE
    XI = 40.0*(HTV/TOP)
    IDX = INT(XI) + 1
    IF (IDX.GE.40) THEN
      CCHHAT = CCH(39)*(40.0 - XI)
    ELSE
      XXI = REAL(IDX)
      CCHHAT = CCH(IDX) + (CCH(IDX-1)-CCH(IDX))*(XXI - XI)
    ENDIF
  ENDIF
  IF (CCHHAT.LT.0.0) CCHHAT = 0.0
  CCHT(I) = CCH_A + CCH_B*CCHHAT
ENDDO
RETURN
END


SUBROUTINE GOMPON
!  GOMPMORT keyword hook (called from vbase/initre.f90 when the GOMPMORT keyword
!  is read). Flags keyword activation; GOMPLOAD (called from MORCON, which runs
!  after the keyword reader) then loads coefficients and sets LGOMP. The coeff
!  file path still comes from FVS_GOMPIT_COEF (a deployment detail); the keyword
!  records the on/off choice reproducibly in the keyfile.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
LGOMPKW = .TRUE.
RETURN
END


BLOCK DATA GOMPBD
!  Initialise the gompit COMMON flags so LGOMP / LGOMPKW / GHAVE are defined
!  before any activation. Required because the morts hook tests
!  IF(LGOMP.AND.GHAVE(ISPC)) and Fortran does not guarantee short-circuit .AND.
IMPLICIT NONE
INCLUDE 'PRGPRM.f90'
INCLUDE 'GOMPMC.f90'
DATA LGOMP /.FALSE./, LGOMPKW /.FALSE./, GHAVE /MAXSP*.FALSE./
DATA NGOMP /0/, GGRP /MAXSP*16/
DATA GDBHMIN /1.0/, GSFLOOR /0.95/
END
