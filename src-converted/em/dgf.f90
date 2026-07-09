SUBROUTINE DGF(DIAM)
IMPLICIT NONE
!----------
! EM $Id$
!----------
!  THIS SUBROUTINE COMPUTES THE VALUE OF DDS (CHANGE IN SQUARED
!  DIAMETER) FOR EACH TREE RECORD, AND LOADS IT INTO THE ARRAY
!  WK2.  DDS IS PREDICTED FROM HABITAT TYPE, LOCATION, SLOPE,
!  ASPECT, ELEVATION, DBH, CROWN RATIO, BASAL AREA IN LARGER TREES,
!  AND CCF.  THE SET OF TREE DIAMETERS TO BE USED IS PASSED AS THE
!  ARGUEMENT DIAM.  THE PROGRAM THUS HAS THE FLEXIBILITY TO
!  PROCESS DIFFERENT CALIBRATION OPTIONS.  THIS ROUTINE IS CALLED
!  BY **DGDRIV** DURING CALIBRATION AND WHILE CYCLING FOR GROWTH
!  PREDICTION.  ENTRY **DGCONS** IS CALLED BY **RCON** TO LOAD SITE
!  DEPENDENT COEFFICIENTS THAT NEED ONLY BE RESOLVED ONCE.
!----------
!OMMONS
!
!
INCLUDE 'PRGPRM.f90'
!
!
INCLUDE 'PLOT.f90'
!
!
INCLUDE 'ARRAYS.f90'
!
!
INCLUDE 'CALCOM.f90'
!
!
INCLUDE 'COEFFS.f90'
!
!
INCLUDE 'CONTRL.f90'
!
!
INCLUDE 'PDEN.f90'
!
!
INCLUDE 'GGCOM.f90'
!
!
INCLUDE 'VARCOM.f90'
INCLUDE 'GREGMC.f90'
!
!
!OMMONS
!----------
!  DIMENSIONS FOR INTERNAL VARIABLES.
!
!     DIAM -- ARRAY LOADED WITH TREE DIAMETERS (PASSED AS AN
!             ARGUEMENT).
!     DGLD -- ARRAY CONTAINING COEFFICIENTS FOR THE LOG(DIAMETER)
!             TERM IN THE DDS MODEL (ONE COEFFICIENT FOR EACH
!             SPECIES).
!     DGCR -- ARRAY CONTAINING THE COEFFICIENTS FOR THE CROWN
!             RATIO TERM IN THE DDS MODEL (ONE COEFFICIENT FOR
!             EACH SPECIES).
!   DGCRSQ -- ARRAY CONTAINING THE COEFFICIENTS FOR THE CROWN
!             RATIO SQUARED TERM IN THE DDS MODEL (ONE
!             COEFFICIENT FOR EACH SPECIES).
!    DGBAL -- ARRAY CONTAINING COEFFICIENTS FOR THE BASAL AREA IN
!             LARGER TREES TERM IN THE DDS MODEL
!             (ONE COEFFICIENT FOR EACH SPECIES).
!   DGDBAL -- ARRAY CONTAINING COEFFICIENTS FOR THE INTERACTION
!             BETWEEN BASAL AREA IN LARGER TREES AND LN(DBH) (ONE
!             COEFFICIENT PER SPECIES).
!    DGCCF -- ARRAY CONTAINING THE COEFFICIENTS FOR THE CROWN
!             COMPETITION FACTOR TERM IN THE DDS MODEL (ONE
!             COEFFICIENT FOR EACH SPECIES, LOADED IN RCON).
!----------
LOGICAL DEBUG
INTEGER K,ISPC,I1,I2,I3,I,IPCCF,IDTYPE,JDTYPE,L
REAL ALCCF,CONSPP,DGLDS,DGBALS,DGCRS,DGCRS2,DGDSQS
REAL DGDBLS,D,ALD,CR,BAL,DUM1,DUM2,DDS,XSLOPE
REAL GGD,GGINC,GGBK,GGDI,GGDS
INTEGER IGGY
REAL DIAM(MAXTRE),DGLD(MAXSP),DGBAL(MAXSP),DGCR(MAXSP), &
      DGCRSQ(MAXSP),DGDBAL(MAXSP),DGHAB(8,MAXSP),DGFOR(6,MAXSP), &
      DGDS(4,MAXSP),DGEL(MAXSP),DGEL2(MAXSP),DGSASP(MAXSP), &
      DGCASP(MAXSP),DGSLOP(MAXSP),DGSLSQ(MAXSP)
INTEGER MAPHAB(117,MAXSP),MAPDSQ(7,MAXSP),MAPLOC(7,MAXSP), &
      MAPCCF(30)
REAL OBSERV(6,MAXSP),DGLCCF(MAXSP),DGPCC1(MAXSP),DGPCCF(MAXSP), &
      DGPCC2(MAXSP),DGCCFA(3)
INTEGER ICLS,LSI,ISIC
REAL YOUNG,OLD,BA100,SI,BARK,DPP,BATEM,DF,DIAGR,ASPDG, &
      XSITE,TMPASP,BRATIO
!----------
!  SPECIES ORDER:
!   1=WB,  2=WL,  3=DF,  4=LM,  5=LL,  6=RM,  7=LP,  8=ES,
!   9=AF, 10=PP, 11=GA, 12=AS, 13=CW, 14=BA, 15=PW, 16=NC,
!  17=PB, 18=OS, 19=OH
!
!  SPECIES EXPANSION
!  LM USES IE LM (ORIGINALLY FROM TT VARIANT)
!  LL USES IE AF (ORIGINALLY FROM NI VARIANT)
!  RM USES IE JU (ORIGINALLY FROM UT VARIANT)
!  AS,PB USE IE AS (ORIGINALLY FROM UT VARIANT)
!  GA,CW,BA,PW,NC,OH USE IE CO (ORIGINALLY FROM CR VARIANT)
!----------
DATA DGLD/ &
     0.80110,  0.80110,  0.82617,  0.213947,   0.86240,        0., &
     0.75719,  0.83323,  0.89571,   0.58932,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,   0.80110, &
          0./
DATA DGCR/ &
     1.02878,  1.02878,  1.41015,  1.523464,   0.52044,        0., &
     1.43611,  1.10040,  1.30934,   1.94874,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,   1.02878, &
          0./
DATA DGCRSQ/ &
    -0.45448, -0.45448, -0.55362,        0.,   0.86236,        0., &
    -0.44926,       0., -0.17348,  -0.88761,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.45448, &
          0./
DATA DGBAL/ &
     0.00064,  0.00064,  0.00216, -0.358634,        0.,        0., &
     0.00556,  0.00443,  0.00350,   0.00710,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,   0.00064, &
          0./
DATA DGDBAL/ &
    -0.00328, -0.00328, -0.00889,        0.,  -0.51270,        0., &
    -0.01282, -0.01281, -0.00764,  -0.02316,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.00328, &
          0./
DATA DGLCCF/ &
    -0.25717, -0.25717, -0.11974,        0.,        0.,        0., &
    -0.12332,       0.,       0.,        0.,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.25717, &
          0./
DATA DGPCC1/ &
         0.,        0., -0.00211,        0.,        0.,        0., &
   -0.00120,  -0.00155, -0.00164,  -0.00297,        0.,        0., &
         0.,        0.,       0.,        0.,        0.,        0., &
         0./
DATA DGPCC2/ &
         0.,        0., -0.00232,        0.,        0.,        0., &
   -0.00145,  -0.00123, -0.00232,  -0.00254,        0.,        0., &
         0.,        0.,       0.,        0.,        0.,        0., &
         0./
DATA DGPCCF/ &
          0.,       0.,  0.000003,       0.,       0.,       0., &
    0.000002, 0.000002, -0.000003, 0.000005,       0.,       0., &
          0.,       0.,        0.,       0.,       0.,       0., &
          0./
!
DATA  OBSERV/ &
     256.,    0.,    0.,    0.,    0.,    0., &
     306.,    0.,    0.,    0.,    0.,    0., &
     697.,    0.,    0.,    0.,    0.,    0., &
      27.,   70.,  123.,  101.,   33.,    0., &
     988.,  177.,   40.,  331.,  332., 2380., &
    1000.,    0.,    0.,    0.,    0.,    0., &
    3208.,    0.,    0.,    0.,    0.,    0., &
       0.,    0.,    0.,    0.,    0.,    0., &
    1198.,    0.,    0.,    0.,    0.,    0., &
    7837.,    0.,    0.,    0.,    0.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0., &
     184.,  429.,  356.,  162.,   74.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0., &
     184.,  429.,  356.,  162.,   74.,    0., &
      22.,    0.,    0.,    0.,    0.,    0., &
    1000.,    0.,    0.,    0.,    0.,    0./
!----------
!  DGHAB IS AN ARRAY THAT CONTAINS HABITAT CLASS INTERCEPTS FOR
!  EACH SPECIES.  MAPHAB IS INDEXED BY JDTYPE TO MAP HABITAT
!  TYPES ONTO HABITAT CLASSES.
!----------
DATA ((MAPHAB(L,K),L=1,117),K=1,3)/ &
   2*4,4*3,16*2,5*4,2*3,2*4,4*3,4,9*3,4,2*3,4*4,2*3,2*2, &
   3*3,2,2*3,2,13*3,3*2,2*3,2*2,3*3,1,8*3,2,5*3,4,3,2,3, &
   2,3,2,3,2*2,5*3,4, &
! SPECIES 2=WL
   2*4,4*3,16*2,5*4,2*3,2*4,4*3,4,9*3,4,2*3,4*4,2*3,2*2, &
   3*3,2,2*3,2,13*3,3*2,2*3,2*2,3*3,1,8*3,2,5*3,4,3,2,3, &
   2,3,2,3,2*2,5*3,4, &
! SPECIES 3=DF
   2*5,7,2*4,2*7,3,3*5,6,12*4,2,4,2*5,2*3,4,2*5,3,7,3*5,7,5, &
   2*4,2*5,2,6,5,5*3,3*7,2*5,4*6,3*5,3,7,2*3,2*6,5*3,5,2*6, &
   4*5,4,3,2*5,3*7,5,2*6,5*5,4,6,1,5,2*3,8*5,6*4/
!
DATA ((MAPHAB(L,K),L=1,117),K=4,6)/ &
   30*1,87*0, &
! SPECIES 5=LL
     12*6, 1,2,3,3,4,3,1, 7*6, 1,5,5,6,87*0, &
! SPECIES 6=RM
   30*1,87*0/
!
DATA ((MAPHAB(L,K),L=1,117),K=7,9)/ &
   1,7,6,1,2,2*6,8,7,2*1,8,3,8,2,1,8,1,3,3*1,2*8,2*3,2,7,3,7, &
   2*1,2,3,1,2*3,2*1,7,3,7,2*1,2,2*3,1,3*2,3,7,6,2*1,2,7,3,2*8, &
   3*7,8,6,8,6,1,7,3,6,2,3,4*2,7,2,3,1,7,3*2,5,1,2,7,8,2,4,2,5, &
   2,1,2*3,3*1,3,1,4*6,2,3*6,5,3,6,3,7, &
! SPECIES 8=ES
   24*4,5,1,4,5,2*6,3,2*4,2*5,4*4,3*3,3*4,2,5,6,3*4, &
   5,5*4,3*5,4,5,2*4,2*2,1,2*4,1,2,4,3*5,3,4,2*2,4,2, &
   2*4,2,3,3*4,3,3*4,3,2*4,3*2,6,5,2*4,3,5,4,3,5*4,6*6, &
! SPECIES 9=AF
   3,27*2,3*5,8*2,13*5,2,3*5,10*4,5*3,1,3,5*2,4,2,2*4, &
   2*3,5*2,2*3,4,2,3,2,3,4,5*3,5,3,7*2,6*5/
!
DATA ((MAPHAB(L,K),L=1,117),K=10,12)/ &
   12*1,3*2,3,3,2,3,4,5,4,1,3*2,6,3,3,10*4,9*5,6*7, &
   3,3,61*7, &
! SPECIES 11=GA
   30*1,87*0, &
! SPECIES 12=AS
   30*1,87*0/
!
DATA ((MAPHAB(L,K),L=1,117),K=13,15)/ &
   30*1,87*0, &
! SPECIES 14=BA
   30*1,87*0, &
! SPECIES 15=PW
   30*1,87*0/
!
DATA ((MAPHAB(L,K),L=1,117),K=16,MAXSP)/ &
   30*1,87*0, &
! SPECIES 17=PB
   30*1,87*0, &
! SPECIES 18=OS
   2*4,4*3,16*2,5*4,2*3,2*4,4*3,4,9*3,4,2*3,4*4,2*3,2*2, &
   3*3,2,2*3,2,13*3,3*2,2*3,2*2,3*3,1,8*3,2,5*3,4,3,2,3, &
   2,3,2,3,2*2,5*3,4, &
! SPECIES 19=OH
   30*1,87*0/
!
DATA ((DGHAB(L,K),L=1,8),K=1,10)/ &
    0.,-0.20842, 0.01545, 0.29742,      0.,      0.,      0.,     0., &
! SPECIES 2=WL
    0.,-0.20842, 0.01545, 0.29742,      0.,      0.,      0.,     0., &
! SPECIES 3=DF
    0.,-0.20240,-0.03260, 0.05559,-0.12847,-0.07100,-0.31860,     0., &
! SPECIES 4=LM
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 5=LL
   -0.96389,-0.72415,-0.57308,-0.82218,-1.24093,-1.10746,  0.,    0., &
! SPECIES 6=RM
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 7=LP
    0.,-0.06686, 0.00921,-0.17011,-0.11806,-0.22504, 0.06528,0.14340, &
! SPECIES 8=ES
    0., 0.23303,-0.01746, 0.11129, 0.19065, 0.40362,      0.,     0., &
! SPECIES 9=AF
    0., 0.12967, 0.33768, 0.25940, 0.48657,      0.,      0.,     0., &
! SPECIES 10=PP
    0.,-0.19356,-0.07999, 0.03500,-0.13880,-1.04929,-0.32378,     0./
!
DATA ((DGHAB(L,K),L=1,8),K=11,MAXSP)/ &
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 12=AS
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 13=CW
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 14=BA
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 15=PW
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 16=NC
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 17=PB
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0., &
! SPECIES 18=OS
    0.,-0.20842, 0.01545, 0.29742,      0.,      0.,      0.,     0., &
! SPECIES 19=OH
    0.,      0.,      0.,      0.,      0.,      0.,      0.,     0./
!----------
!  FOR SUBALPINE LARCH FROM THE IE VARIANT,
!  DGCCFA CONTAINS COEFFICIENTS FOR THE CCF TERM BY SPECIES BY
!  HABITAT CLASS.  MAPCCF IS INDEXED BY ITYPE TO MAP HABITAT TYPES
!  ONTO HABITAT CLASSES.
!  LIMBER PINE USES THE SAME COEFFICIENT FOR ALL HABITAT TYPES AND IS
!  LOADED IN ENTRY DGCONS
!----------
DATA MAPCCF/12*3, 7*1, 3,1,1,2,2,1,3,1,1,3,3/
DATA DGCCFA/ -0.01598, -0.04477, -0.07392/
!----------
!  DGFOR CONTAINS LOCATION CLASS CONSTANTS FOR EACH SPECIES.
!  MAPLOC IS AN ARRAY WHICH MAPS FOREST ONTO A LOCATION CLASS.
!----------
!
DATA ((MAPLOC(L,K),L=1,7),K=1,10)/ &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,2,3,4,3,1,5, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,2,3,2,4,5,6, &
    1,2,3,2,3,4,1, &
    1,2,3,4,5,6,1, &
    1,2,3,3,2,2,1/
!
DATA ((MAPLOC(L,K),L=1,7),K=11,MAXSP)/ &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1/
!---------
! WB,WL AND OT WERE SHOWING TOO SLOW DIAMETER GROWTH.
! CHANGED LOCATION CONSTANT FOR WB,WL,OT FROM -0.15675
! TO 1.5675. 07-11-01
!---------
!
DATA ((DGFOR(L,K),L=1,6),K=1,10)/ &
     1.5675,      0.,      0.,      0.,      0.,      0., &
     1.5675,      0.,      0.,      0.,      0.,      0., &
    1.10349, 1.55392, 1.04495, 1.27679,      0.,      0., &
   1.568742,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
    1.75497, 1.89439, 1.63227, 1.65474, 1.51043,      0., &
   -2.45844,-2.14138,-2.55288,-2.42650,      0.,      0., &
    0.46110,-0.16517, 0.35319, 0.64719, 0.40140, 0.33757, &
    3.57069, 3.45044, 3.67551,      0.,      0.,      0./
!
DATA ((DGFOR(L,K),L=1,6),K=11,MAXSP)/ &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0., &
     1.5675,      0.,      0.,      0.,      0.,      0., &
         0.,      0.,      0.,      0.,      0.,      0./
!----------
!  DGDS CONTAINS COEFFICIENTS FOR THE DIAMETER SQUARED TERMS
!  IN THE DIAMETER INCREMENT MODELS; ARRAYED BY FOREST BY
!  SPECIES.  MAPDSQ IS AN ARRAY WHICH MAPS FOREST ONTO A DBH**2
!  COEFFICIENT.
!----------
!
DATA ((MAPDSQ(L,K),L=1,7),K=1,10)/ &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,3,2,1,3,2,3, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,2, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,2,1,1,1,1,1/
!
DATA ((MAPDSQ(L,K),L=1,7),K=11,MAXSP)/ &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1, &
    1,1,1,1,1,1,1/
!
DATA ((DGDS(L,K),L=1,4),K=1,10)/ &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
     -0.000251, -0.000412,            0.,       0., &
    -0.0006538,        0.,            0.,       0., &
     -0.000283,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
     -0.000402,        0.,            0.,       0., &
     -0.000034,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
     -0.000248, -0.000991,            0.,       0./
!
DATA ((DGDS(L,K),L=1,4),K=11,MAXSP)/ &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0., &
            0.,        0.,            0.,       0./
!----------
!  DGEL CONTAINS THE COEFFICIENTS FOR THE ELEVATION TERM IN THE
!  DIAMETER GROWTH EQUATION.  DGEL2 CONTAINS THE COEFFICIENTS FOR
!  THE ELEVATION SQUARED TERM IN THE DIAMETER GROWTH EQUATION.
!  DGSASP CONTAINS THE COEFFICIENTS FOR THE SIN(ASPECT)*SLOPE
!  TERM IN THE DIAMETER GROWTH EQUATION.  DGCASP CONTAINS THE
!  COEFFICIENTS FOR THE COS(ASPECT)*SLOPE TERM IN THE DIAMETER
!  GROWTH EQUATION.  DGSLOP CONTAINS THE COEFFICIENTS FOR THE
!  SLOPE TERM IN THE DIAMETER GROWTH EQUATION.  DGSLSQ CONTAINS
!  COEFFICIENTS FOR THE (SLOPE)**2 TERM IN THE DIAMETER GROWTH MODELS.
!  ALL OF THESE ARRAYS ARE SUBSCRIPTED BY SPECIES.
!----------
DATA DGCASP/ &
     0.00270,  0.00270, -0.04290, -0.609774,  -0.12473,        0., &
    -0.02355, -0.15989, -0.20925,   0.06497,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,   0.00270, &
          0./
DATA DGSASP/ &
    -0.01606, -0.01606,  0.06124,  -0.01752,  -0.06862,        0., &
    -0.00393,  0.03982,  0.01463,   0.20328,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.01606, &
          0./
DATA DGSLOP/ &
    -0.20011, -0.20011, -0.39500,  -2.05706,   0.30070,        0., &
    -0.21964,  0.06898, -0.03350,  -0.95238,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.20011, &
          0./
DATA DGSLSQ/ &
          0.,       0.,       0.,  2.113263,  -0.62224,        0., &
          0., -0.34251, -0.11164,   0.61813,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,        0., &
          0./
DATA DGEL/ &
    -0.00565, -0.00565, -0.00196,        0.,   0.06313,        0., &
    -0.01234,  0.08677, -0.00885,  -0.07453,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,  -0.00565, &
          0./
DATA DGEL2/ &
          0.,       0.,       0.,        0., -0.000676,        0., &
          0., -0.00072, -0.00002,   0.00060,        0.,        0., &
          0.,       0.,       0.,        0.,        0.,        0., &
          0./
!-----------
!  CHECK FOR DEBUG.
!-----------
CALL DBCHK (DEBUG,'DGF',3,ICYC)
!
IF(DEBUG) WRITE(JOSTND,3)ICYC
3 FORMAT(' ENTERING SUBROUTINE DGF  CYCLE =',I5)
!----------
!  DEBUG OUTPUT: MODEL COEFFICIENTS.
!----------
IF(DEBUG) &
    WRITE(JOSTND,9000) DGCON,DGDSQ,DGLD,DGCR,DGCRSQ,DGCCF,DGBAL
9000 FORMAT(/,' DGCON',/12(1X,F9.5),/7(1X,F9.5),/,' DGDSQ',/12(1X,F9.5) &
    ,/7(1X,F9.5),/,' DGLD',/12(1X,F9.5),/7(1X,F9.5),/,' DGCR', &
    /12(1X,F9.5),/7(1X,F9.5),/,' DGCRSQ',/12(1X,F9.5),/7(1X,F9.5), &
    /,' DGCCF',/12(1X,F9.5),/7(1X,F9.5),/,' DGBAL',/12(1X,F9.5), &
    /7(1X,F9.5))
!----------
! CALL BADIST AND COMPUTE THE AGE RANGE
!----------
AGERNG= 1000.
YOUNG= 1000.
OLD= 0.0
CALL BADIST(DEBUG)
!
DO I= 1,ITRN
  IF(ABIRTH(I).LE.1.)GO TO 5
  IF(HT(I) .LE. 4.5) GO TO 5
  IF(ABIRTH(I) .LT. YOUNG) YOUNG = ABIRTH(I)
  IF(ABIRTH(I) .GT. OLD) OLD = ABIRTH(I)
5   CONTINUE
  AGERNG = OLD - YOUNG
  AGERNG = ABS(AGERNG)
ENDDO
!
CALL SDICAL(0,SDIMAX)
!----------
!  SCALE BASAL AREA.
!----------
BA100 = BA/100.
!----------
!  BEGIN SPECIES LOOP.  ASSIGN VARIABLES WHICH ARE SPECIES DEPENDENT
!----------
ALCCF=0.0
IF(RELDEN .GT. 0.0)ALCCF=ALOG(RELDEN)
CALL GREGLOADDG
DO 20 ISPC=1,MAXSP
I1=ISCT(ISPC,1)
IF(I1.EQ.0) GO TO 20
I2=ISCT(ISPC,2)
CONSPP= DGCON(ISPC) + COR(ISPC) + 0.01*DGCCF(ISPC)*RELDEN
!
!  SPECIES EXPANSION-CONSPP IS DIFFERENT FOR LM
!
IF(ISPC.EQ.4)CONSPP= CONSPP + .01*(-.199592)*RELDEN
!
DGLDS= DGLD(ISPC)
DGBALS = DGBAL(ISPC)
DGCRS= DGCR(ISPC)
DGCRS2=DGCRSQ(ISPC)
DGDSQS=DGDSQ(ISPC)
DGDBLS=DGDBAL(ISPC)
!
!  SPECIES EXPANSION-SITE INDEX ASSIGNED IS DEPENDENT UPON HABITAT TYPE
!  IN SUBROUTINE **SITSET**
!
SI = SITEAR(ISPC)
!----------
!  BEGIN TREE LOOP WITHIN SPECIES ISPC.
!----------
DO 10 I3=I1,I2
I=IND1(I3)
D=DIAM(I)
IF (D.LE.0.0) GOTO 10
ALD=ALOG(D)
BARK=BRATIO(ISPC,D,HT(I))
CR=ICR(I)*0.01
BAL = (1.0 - (PCT(I)/100.)) * BA
IPCCF=ITRE(I)
!
!  RM JUNIPER LOGIC FROM UT
!
IF(ISPC.EQ.6)THEN
  DPP = D
  IF(DPP .LT. 1.0) DPP = 1.0
  BATEM = BA
  IF(BATEM .LT.  1.0) BATEM = 1.0
  DF = 0.25897 + 1.03129 * DPP - 0.0002025464 * BATEM &
            + 0.00177 * SI
  IF((DF-DPP) .GT. 1.0) DF = DPP + 1.0
  IF(DF .LT. DPP) DF = DPP
  DIAGR = (DF - DPP) * BARK
  IF(DIAGR .LE. 0.) THEN
    DDS=-9.21
  ELSE
    DDS = ALOG( (DIAGR * (2.0 * DPP * BARK + DIAGR)) ) + CONSPP
  ENDIF
!
!  ASPEN LOGIC FROM UTAH; ALSO USED FOR PAPER BIRCH
!
ELSEIF(ISPC.EQ.12 .OR. ISPC.EQ.17)THEN
  CR=ICR(I)
  CALL DGFASP(D,ASPDG,CR,BARK,SI,DEBUG)
  DDS = ASPDG + ALOG(COR2(ISPC)) + COR(ISPC)
!
!  CO LOGIC FROM CR
!
ELSEIF(ISPC.EQ.11 .OR. (ISPC.GE.13 .AND. ISPC.LE.16) .OR. &
          ISPC.EQ.19) THEN
  ICLS = IFIX(D + 1.0)
  IF(ICLS .GT. 41) ICLS = 41
  DPP = D
  IF(DPP .LT. 1.0) DPP = 1.0
  DF = 0.24506 + 1.01291 * DPP - 0.00084659 * BA &
            + 0.00631*SI
  IF(DF .GT. 36.) DF = 36.
  IF(DF .LT. DPP) DF = DPP
  DIAGR = (DF - DPP) * BARK
  IF(DIAGR .LE. 0.) THEN
    DDS=-9.21
  ELSE
    DDS = ALOG( (DIAGR * (2.0 * DPP * BARK + DIAGR)) )
    IF(DDS .LT. -9.21) DDS=-9.21
  ENDIF
  DDS = DDS+COR(ISPC)+DGCON(ISPC)
!
! ORIGINAL NI SECTION
!
ELSEIF(ISPC.EQ.4 .OR. ISPC.EQ.5)THEN
  BAL = (1.0-(PCT(I)/100.))*BA100
  DDS=CONSPP + DGLDS*ALD + DGBALS*BAL + CR*(DGCRS+CR*DGCRS2) &
             +DGDSQS*D*D  + DGDBLS*BAL/(ALOG(D+1.0))
ELSE
!----------
! ORIGINAL EM SECTION
! DUM1 AND DUM2 ARE RELDEN DUMMYS.  THEY ARE 0 OR 1 DEPENDING ON WHETHER
! OR NOT THE STAND IS MANAGED.
!----------
  DUM1 = 0.0
  DUM2 = 1.0
  IF(MANAGD .EQ. 1) THEN
    DUM1 = 1.0
    DUM2 = 0.0
  ENDIF
  DDS=CONSPP + DGLDS*ALD + DGBALS*BAL + CR*(DGCRS+CR*DGCRS2) &
      + DGDSQS*D*D  + DGDBLS*BAL/(ALOG(D+1.0)) &
      + DGPCCF(ISPC)*(RELDEN**2) &
      + DGLCCF(ISPC) * ALCCF + DGPCC1(ISPC)*DUM1*PCCF(IPCCF) &
      + DGPCC2(ISPC)*DUM2*PCCF(IPCCF)
!
  IF(DEBUG) WRITE(JOSTND,8000) DGCON(ISPC),COR(ISPC), &
     DGCCF(ISPC),RELDEN,DGLCCF(ISPC),RMAI,CONSPP,DGLDS,ALD, &
     DGBALS,BAL,CR,DGCRS,DGCRS2,DGDSQS,D,DGDBLS, &
     DGPCCF(ISPC),PCCF(IPCCF),DDS,DGLCCF(ISPC),ALCCF
8000   FORMAT(/,' IN DGF AT LINE 178',9F12.5,/,1H ,8F12.5,/, &
     1H ,5F12.5,//)
ENDIF
!---------
IF(DDS.LT.-9.21) DDS=-9.21
WK2(I)=DDS
      IF (LGREGDG .AND. GHAVE_DG(ISPC)) THEN
        GGD = DIAM(I)
        DO IGGY = 1, 10
          CALL GREGDGV(ISPC, GGD, FLOAT(ICR(I))/100.0, HT(I), BAL, GELEV, GEMT, GGINC)
          GGD = GGD + GGINC
        END DO
        GGBK = BRATIO(ISPC, GGD, HT(I))
        GGDI = (GGD - DIAM(I)) * GGBK
        GGDS = GGDI * (2.0*DIAM(I)*GGBK + GGDI)
        IF (GGDS .LT. 1.0E-6) GGDS = 1.0E-6
        WK2(I) = ALOG(GGDS)
      ENDIF
!----------
!  END OF TREE LOOP.  PRINT DEBUG INFO IF DESIRED.
!----------
IF(DEBUG)THEN
   WRITE(JOSTND,9001) I,ISPC,D,BAL,CR,RELDEN,BA,DDS
9001    FORMAT(' IN DGF, I=',I4,',  ISPC=',I3,',  DBH=',F7.2, &
         ',  BAL=',F7.2,',  CR=',F7.4/ &
         '       RELDEN=',F9.3,',  BA=',F9.3,',   LN(DDS)=',F7.4)
ENDIF
10 CONTINUE
!----------
!  END OF SPECIES LOOP.
!----------
20 CONTINUE
IF(DEBUG)WRITE(JOSTND,7000)ICYC
7000 FORMAT(' LEAVING SUBROUTING DGF  CYCLE =',I5)
RETURN
!
!
!
ENTRY DGCONS(IDTYPE)
!----------
!  ENTRY POINT FOR LOADING COEFFICIENTS OF THE DIAMETER INCREMENT
!  MODEL THAT ARE SITE SPECIFIC AND NEED ONLY BE RESOLVED ONCE.
!  IDTYPE IS A HABITAT TYPE INDEX THAT IS COMPUTED IN **RCON**.
!  ASPECT IS STAND ASPECT.  OBSERV CONTAINS THE NUMBER OF
!  OBSERVATIONS BY HABITAT CLASS BY SPECIES FOR THE UNDERLYING
!  MODEL (THIS DATA IS ACTUALLY USED BY **DGDRIV** FOR CALIBRATION).
!----------
JDTYPE=IDTYPE
IF(IDTYPE .GT. 117) JDTYPE=30
!----------
!  ENTER LOOP TO LOAD SPECIES DEPENDENT VECTORS.
!----------
DO 30 ISPC=1,MAXSP
!------------
!  PUT SITE INTO ONE OF 5 CLASSES
!----------
XSITE = SITEAR(ISPC)
LSI=INT(XSITE/10.)
IF(LSI.LT.2)ISIC=1
IF(LSI.GE.5)ISIC=5
IF(LSI.GE.3 .AND. LSI.LT.4)ISIC=3
IF(LSI.GE.4 .AND. LSI.LT.5)ISIC=4
IF(LSI.GE.2 .AND. LSI.LT.3)ISIC=2
!
IF(ISPC.LE.3 .OR. (ISPC.GE.7 .AND. ISPC.LE.10) .OR. &
   ISPC.EQ.18)THEN
  ISPHAB=MAPHAB(JDTYPE,ISPC)
ELSE
  ISPHAB=MAPHAB(ITYPE,ISPC)
ENDIF
ISPFOR=MAPLOC(IFOR,ISPC)
ISPDSQ=MAPDSQ(IFOR,ISPC)
ISPCCF=MAPCCF(ITYPE)
TMPASP = ASPECT
IF(ISPC.EQ.4 .OR. ISPC.EQ.6 .OR. ISPC.EQ.12 &
      .OR.ISPC.EQ.17)TMPASP=TMPASP-.7854
XSLOPE=SLOPE
IF(ISPC.LE.3 .OR. (ISPC.GE.7 .AND. ISPC.LE.10) .OR. &
      ISPC.EQ.18) XSLOPE=SLOPE / 10.0
!
! NULLIFY THE SLOPE COEF FOR OS BECAUSE OF POOR EXTRAPOLATION
!
IF(ISPC .EQ. 18)XSLOPE=0.0
DGCON(ISPC)= DGHAB(ISPHAB,ISPC) &
                    + DGFOR(ISPFOR,ISPC) &
                    + DGEL(ISPC) * ELEV &
                    + DGEL2(ISPC) * ELEV * ELEV &
                    +(DGSASP(ISPC) * SIN(TMPASP) &
                    + DGCASP(ISPC) * COS(TMPASP) &
                    + DGSLOP(ISPC)) * XSLOPE &
                    + DGSLSQ(ISPC) * XSLOPE * XSLOPE
DGDSQ(ISPC)=DGDS(ISPDSQ,ISPC)
DGCCF(ISPC)=0.0
IF(ISPC.EQ.4)THEN
    DGCCF(ISPC)= -0.199592
    ATTEN(ISPC) = OBSERV(ISIC,ISPC)
ELSEIF(ISPC.EQ.5)THEN
    DGCCF(ISPC)=DGCCFA(ISPCCF)
    ATTEN(ISPC)=OBSERV(ISPHAB,ISPC)
ELSEIF(ISPC.EQ.12 .OR. ISPC.EQ.17)THEN
  ATTEN(ISPC) = OBSERV(ISIC,ISPC)
ELSE
  ATTEN(ISPC)=OBSERV(1,ISPC)
ENDIF
IF(ISPC.EQ.4)THEN
   DGCON(ISPC)=DGCON(ISPC)+.001766*XSITE
ELSEIF(ISPC.EQ.12 .OR. ISPC.EQ.17)THEN
   DGCON(ISPC)=DGCON(ISPC)+.006460*XSITE
ENDIF
SMCON(ISPC)=0.
!----------
!  IF READCORD OR REUSCORD WAS SPECIFIED (LDCOR2 IS TRUE) ADD
!  LN(COR2) TO THE BAI MODEL CONSTANT TERM (DGCON).  COR2 IS
!  INITIALIZED TO 1.0 IN BLKDATA.
!----------
IF (LDCOR2.AND.COR2(ISPC).GT.0.0) DGCON(ISPC)=DGCON(ISPC) &
      + ALOG(COR2(ISPC))
30 CONTINUE
RETURN
END
