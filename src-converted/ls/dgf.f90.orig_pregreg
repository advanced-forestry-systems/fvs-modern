SUBROUTINE DGF(DIAM)
IMPLICIT NONE
!----------
! LS $Id$
!----------
!  THIS SUBROUTINE COMPUTES THE VALUE OF DDS (CHANGE IN SQUARED
!  DIAMETER) FOR EACH TREE RECORD, AND LOADS IT INTO THE ARRAY
!  WK2.  DIAMETER GROWTH IS PREDICTED FROM DBH, SITE INDEX,
!  CROWN RATIO, AND BASAL AREA PERCENTILES AND QMD FOR TREES
!  GREATER THAN OR EQUAL TO 5 INCHED DBH. THE
!  SET OF TREE DIAMETERS TO BE USED IS PASSED AS THE ARGUMENT TO
!  DIAM.  THE PROGRAM THUS HAS THE FLEXIBILITY TO PROCESS DIFFERENT
!  CALIBRATION OPTIONS.  THIS ROUTINE IS CALLED BY **DGDRIV** DURING
!  CALIBRATION AND WHILE CYCLING FOR GROWTH PREDICTION.  ENTRY
!  **DGCONS** IS CALLED BY **RCON** TO LOAD SITE DEPENDENT COEFFICIENTS
!  THAT NEED ONLY BE RESOLVED ONCE.
!
!  DIAMETER GROWTH EQUATIONS WERE DEVELOPED THROUGH COST SHARE AGREEMENT NO.
!  10-CS-11132425-258 WITH MICHIGAN TECHNOLOGICAL UNIVERSITY, OCTOBER 2012.
!----------
!OMMONS
!
!
INCLUDE 'PRGPRM.f90'
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
INCLUDE 'OUTCOM.f90'
!
!
INCLUDE 'PDEN.f90'
!
!
INCLUDE 'PLOT.f90'
!
!
!OMMONS
!-----------
!  DEFINITIONS OF LOCAL VARIABLES.
!
!    BAGE5 -- STAND BASAL AREA IN TREES 5.0" DBH AND LARGER.
!      BAL -- STAND BASAL AREA IN TREES WITH DIAMETER >= THE DIAMETER
!             OF THE SUBJECT TREE.
!     BALC -- ARRAY CONTAINING COEFFICIENTS FOR THE BASAL AREA IN
!             LARGER TREES TERM IN THE DDS MODEL.
!     BARK -- BARK RATIO
!   BRATIO -- BARK RATIO FUNCTION CONTAINED IN BRATIO.F
!   CONSPP -- CONSTANT TERM IN THE DIAMETER GROWTH EQUATION.
!       CR -- PERCENT CROWN RATIO FOR A TREE
!    CRSQC -- ARRAY CONTAINING THE COEFFICIENTS FOR THE PERCENT CROWN
!             RATIO SQUARED TERM IN THE DDS MODEL.
!    CRWNC -- ARRAY CONTAINING THE COEFFICIENTS FOR THE PERCENT CROWN
!             RATIO TERM IN THE DDS MODEL.
!        D -- DIAMETER AT BREAST HEIGHT FOR A TREE
!     DBHC -- ARRAY CONTAINING COEFFICIENTS FOR THE DIAMETER TERM IN
!             THE DDS MODEL
!    DBH2C -- ARRAY CONTAINING COEFFICIENTS FOR THE DIAMETER SQUARED
!             TERM IN THE DDS MODEL
!      DDS -- PREDICTED INSIDE-BARK BASAL AREA INCREMENT FOR A TREE
!    DEBUG -- LOGICAL VARIABLE USED TO CONTROL DEBUG OUTPUT
!   DIAGRI -- PREDICTED INSIDE-BARK DIAMETER GROWTH FOR A TREE
!   DIAGRO -- PREDICTED OUTSIDE-BARK DIAMETER GROWTH FOR A TREE
!     DIAM -- ARRAY LOADED WITH TREE DIAMETERS
! D1 - D11 -- EQUATION COEFFICIENTS FOR A DEBUG WRITE STATEMENT
!     I,IK -- ARRAY INDICIES
!   INTERC -- DDS REGRESSION EQUATION INTERCEPT
!     ISPC -- FVS SPECIES SEQUENCE NUMBER
! I1,I2,I3 -- ARRAY INDICIES
!   OBSERV -- SAMPLE SIZE USED TO FIT THE EQUATION
!        P -- TREES-PER-ACRE FOR A TREE RECORD
!   QMDGE5 -- QUADRATIC MEAN DIAMETER FOR TREES 5.0" DBH AND LARGER
!    RDBHC -- ARRAY CONTAINING COEFFICIENTS FOR DIAM/QMDGE5
!             TERM IN THE DDS MODEL.
!  RDBHSQC -- ARRAY CONTAINING COEFFICIENTS FOR THE DIAM*DIAM/QMDGE5
!             TERM IN THE DDS MODEL.
!   RELDBH -- RELATIVE TREE DIAMETER (D/QMDGE5)
! RELDBHSQ -- RELATIVE TREE DIAMTER SQUARED ((D**2)/QMDGE5)
!     SBAC -- ARRAY CONTAINING THE COEFFICIENTS FOR STAND BASAL AREA
!             IN TREES 5.0" DBH AND LARGER TERM IN THE DDS MODEL.
!   SDQGE5 -- SUM OF TREE DIAMETERS SQUARED FOR THE STAND
!    SITEC -- ARRAY CONTAINING COEFFICIENTS FOR SITE INDEX IN THE DDS MODEL.
!        T -- TOTAL STAND TREES-PER-ACRE
!    VDBHC -- ARRAY CONTAINING COEFFICIENTS FOR THE INVERSE OF DIAMETER
!             TERM IN THE DDS MODEL.
!----------
!  SPECIES AND COEFFICIENT ORDER
!----------
!  NOTE -- commercial hardwood (44=CH) uses black walnut (46=WN)
!          as surrogate for all coefficients.
!
!      1=JP   2=SC    3=RN    4=RP    5=WP
!      6=WS   7=NS    8=BF    9=BS   10=TA
!     11=WC  12=EH   13=OS   14=RC   15=BA
!     16=GA  17=EC   18=SV   19=RM   20=BC
!     21=AE  22=RL   23=RE   24=YB   25=BW
!     26=SM  27=BM   28=AB   29=WA   30=WO
!     31=SW  32=BR   33=CK   34=RO   35=BO
!     36=NP  37=BH   38=PH   39=SH   40=BT
!     41=QA  42=BP   43=PB   44=CH   45=BN
!     46=WN  47=HH   48=BK   49=NC   50=BE
!     51=ST  52=MM   53=AH   54=AC   55=HK
!     56=DW  57=HT   58=AP   59=BG   60=SY
!     61=PR  62=CC   63=PL   64=WI   65=BL
!     66=DM  67=SS   68=MA
!-----------
!
LOGICAL DEBUG
!
INTEGER I,IK,ISPC,I1,I2,I3
!
REAL BAGE5,BAL,BARK,BRATIO,CONSPP,CR,D,DDS,DIAGRO,DIAGRI
REAL D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11
REAL P,QMDGE5,RELDBH,RELDBHSQ,SDQGE5,T
!
REAL BALC(MAXSP),CRSQC(MAXSP),CRWNC(MAXSP),DBHC(MAXSP), &
        DBH2C(MAXSP),DIAM(MAXTRE),INTERC(MAXSP), &
        OBSERV(MAXSP),RDBHC(MAXSP), &
        RDBHSQC(MAXSP),SBAC(MAXSP),SITEC(MAXSP), &
        VDBHC(MAXSP)
!----------
!  DATA STATEMENTS:
!
!  COEFFICIENTS FOR DDS MODEL
!----------
DATA INTERC/ &
       1.395,     1.395,  0.21748,   0.21748,   4.1821, &
      3.9303,    3.9303,    3.458,    3.4206,   -0.461, &
    -0.74827,    3.2309,  0.21748,    3.0982,   2.5313, &
      2.5313,    1.8151,   2.2302,   -0.3258,  -0.3258, &
      3.2435,    3.2435,   3.2435,      2.54,   3.6912, &
      3.2946,    3.2946,   2.7284,   2.47033,   1.9019, &
      1.9019,    1.9019,   1.9019,    0.4452,  0.49086, &
     0.49086,  0.058969, 0.058969,  0.058969,   3.4483, &
      3.3904,    3.3904,   2.6884,    2.9591,   2.9591, &
      2.9591,    5.3879,    3.022,    3.1876,   3.1876, &
      3.1876,    3.1876,   3.0308,    3.0308,   1.2352, &
      3.0308,    3.0308,   3.0308,  -0.28233,   3.8669, &
      3.0308,    3.0308,   3.0308,    3.8669,   3.8669, &
      3.8669,    3.3841,   3.0308/
!
DATA VDBHC/ &
    -3.9398,  -3.9398,      0.0,      0.0,  -13.911, &
    -9.3995,  -9.3995,  -10.103,  -9.4927,      0.0, &
        0.0,  -11.816,      0.0,  -8.5844,  -10.146, &
    -10.146,      0.0,  -6.9962,      0.0,      0.0, &
    -10.737,  -10.737,  -10.737,   -7.304,  -11.035, &
    -12.691,  -12.691,  -12.059,  -9.9428,  -5.6950, &
    -5.6950,  -5.6950,  -5.6950,      0.0,      0.0, &
        0.0,      0.0,      0.0,      0.0,  -8.1878, &
    -7.1815,  -7.1815,  -7.3181,  -9.5199,  -9.5199, &
    -9.5199,  -13.259,  -13.627,  -10.287,  -10.287, &
    -10.287,  -10.287,  -8.1059,  -8.1059,  -5.2144, &
    -8.1059,  -8.1059,  -8.1059,      0.0,  -8.8704, &
    -8.1059,  -8.1059,  -8.1059,  -8.8704,  -8.8704, &
    -8.8704,  -10.282,  -8.1059/
!
DATA DBHC/ &
         0.0,      0.0,  0.30253,  0.30253,      0.0, &
         0.0,      0.0,      0.0,      0.0,  0.42230, &
     0.33025,      0.0,  0.30253,      0.0,      0.0, &
         0.0, 0.091526,  0.11324,  0.33549,  0.33549, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0, 0.096247, &
    0.096247, 0.096247, 0.096247,  0.22472,  0.23430, &
     0.23430,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,  0.18300, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0/
!
DATA DBH2C/ &
     0.01213,  0.01213, -0.00885, -0.00885,       0.0, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0, -0.00885,      0.0,       0.0, &
         0.0,      0.0, -0.00186, -0.00703,  -0.00703, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0,      0.0, -0.00397,  -0.00371, &
    -0.00371,  0.01986,  0.01986,  0.01986,       0.0, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0,      0.0,      0.0,  -0.00314, &
         0.0,      0.0,      0.0,  0.01581,       0.0, &
         0.0,      0.0,      0.0,      0.0,       0.0, &
         0.0,      0.0,      0.0/
!
DATA RDBHC/ &
     1.37270,  1.37270,      0.0,      0.0, -0.33396, &
    -0.58596, -0.58596, -0.63924, -1.43830,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,  1.02090,      0.0, -0.53201, -0.53201, &
         0.0,      0.0,      0.0,      0.0, -0.34599, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,  4.09180,  4.09180,  4.09180,      0.0, &
    -0.16910, -0.16910,      0.0,      0.0,      0.0, &
         0.0, -2.44892,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,  3.37463,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0/
!
DATA RDBHSQC/ &
    -0.17735, -0.17735,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0, -0.11080, &
    -0.08958,      0.0,      0.0,      0.0,      0.0, &
         0.0, -0.03235,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0, -0.01611, &
    -0.01611, -0.01611, -0.01611,      0.0,      0.0, &
         0.0, -0.33415, -0.33415, -0.33415,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0, -0.26631,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0/
!
DATA SBAC/ &
         0.0,      0.0, -0.00114,  -0.00114,      0.0, &
         0.0,      0.0, -0.00192, 0.0049528,      0.0, &
    -0.00118, -0.00166, -0.00114,  -0.00386, -0.00095, &
    -0.00095,      0.0,      0.0,  -0.00203, -0.00203, &
         0.0,      0.0,      0.0,       0.0,      0.0, &
    -0.00210, -0.00210,      0.0,  -0.00228,      0.0, &
         0.0,      0.0,      0.0,       0.0, -0.00186, &
    -0.00186,      0.0,      0.0,       0.0,      0.0, &
         0.0,      0.0,      0.0,       0.0,      0.0, &
         0.0,      0.0, -0.00288,       0.0,      0.0, &
         0.0,      0.0, -0.00179,  -0.00179, -0.00439, &
    -0.00179, -0.00179, -0.00179,       0.0,      0.0, &
    -0.00179, -0.00179, -0.00179,       0.0,      0.0, &
         0.0,      0.0, -0.00179/
!
DATA BALC/ &
    -0.00428, -0.00428, -0.00306,   -0.00306, -0.00412, &
    -0.00375, -0.00375,      0.0, -0.0067966,      0.0, &
         0.0,      0.0, -0.00306,        0.0,      0.0, &
         0.0, -0.00428, -0.00171,        0.0,      0.0, &
         0.0,      0.0,      0.0,   -0.00223, -0.00259, &
         0.0,      0.0,      0.0,        0.0, -0.00367, &
    -0.00367, -0.00367, -0.00367,   -0.00148,      0.0, &
         0.0,  0.00184,  0.00184,    0.00184, -0.00134, &
    -0.00311, -0.00311,      0.0,   -0.00169, -0.00169, &
    -0.00169, -0.00353,      0.0,        0.0,      0.0, &
         0.0,      0.0,      0.0,        0.0,  0.00331, &
         0.0,      0.0,      0.0,        0.0, -0.00535, &
         0.0,      0.0,      0.0,   -0.00535, -0.00535, &
    -0.00535, -0.00170,      0.0/
!
DATA CRWNC/ &
    0.02997, 0.02997, 0.01947,  0.01947, 0.03382, &
    0.01635, 0.01635, 0.02677, 0.018393,     0.0, &
    0.02758, 0.01535, 0.01947,  0.01628, 0.01677, &
    0.01677,     0.0, 0.03114,  0.03314, 0.03314, &
    0.05258, 0.05258, 0.05258,  0.01369, 0.01932, &
    0.01787, 0.01787, 0.01843,  0.04226, 0.03021, &
    0.03021, 0.03021, 0.03021,      0.0, 0.00770, &
    0.00770,     0.0,     0.0,      0.0, 0.02367, &
    0.01543, 0.01543,     0.0,  0.03012, 0.03012, &
    0.03012,     0.0, 0.04569,  0.04919, 0.04919, &
    0.04919, 0.04919, 0.00489,  0.00489, 0.02814, &
    0.00489, 0.00489, 0.00489,  0.01650, 0.02034, &
    0.00489, 0.00489, 0.00489,  0.02034, 0.02034, &
    0.02034, 0.03201, 0.00489/
!
DATA CRSQC/ &
    -0.00015, -0.00015,      0.0,      0.0, -0.00012, &
         0.0,      0.0, -0.00009,      0.0,  0.00012, &
    -0.00015,      0.0,      0.0,      0.0,      0.0, &
         0.0,  0.00031, -0.00014, -0.00016, -0.00016, &
    -0.00038, -0.00038, -0.00038,      0.0,      0.0, &
         0.0,      0.0,      0.0, -0.00029, -0.00022, &
    -0.00022, -0.00022, -0.00022,  0.00009,      0.0, &
         0.0,  0.00008,  0.00008,  0.00008,      0.0, &
         0.0,      0.0,  0.00020,      0.0,      0.0, &
         0.0,      0.0, -0.00037, -0.00045, -0.00045, &
    -0.00045, -0.00045,      0.0,      0.0, -0.00011, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0,      0.0,      0.0,      0.0,      0.0, &
         0.0, -0.00026,      0.0/
!
DATA SITEC/ &
    0.00385,  0.00385,  0.00581,   0.00581,       0.0, &
        0.0,      0.0,  0.00317, 0.0074198,       0.0, &
    0.00537,  0.00676,  0.00581,       0.0,   0.01285, &
    0.01285,      0.0,      0.0,   0.00520,   0.00520, &
        0.0,      0.0,      0.0,   0.00647,       0.0, &
    0.00496,  0.00496,  0.00856,   0.01152,       0.0, &
        0.0,      0.0,      0.0,   0.01305,   0.00724, &
    0.00724,      0.0,      0.0,       0.0,       0.0, &
    0.00619,  0.00619,  0.00350,   0.00356,   0.00356, &
    0.00356,      0.0,  0.00766,       0.0,       0.0, &
        0.0,      0.0,      0.0,       0.0,   0.00447, &
        0.0,      0.0,      0.0,       0.0,       0.0, &
        0.0,      0.0,      0.0,       0.0,       0.0, &
        0.0,      0.0,      0.0/
!
DATA OBSERV/ &
    1337.0, 1337.0, 3409.0, 3409.0, 1343.0, &
     986.0,  986.0, 2498.0, 2365.0, 1524.0, &
    5618.0,  756.0, 3409.0, 4111.0, 3285.0, &
    3285.0,  416.0, 3300.0, 5133.0, 5133.0, &
     604.0,  604.0,  604.0,  563.0, 1934.0, &
    4470.0, 4470.0,  150.0,  303.0, 1848.0, &
    1848.0, 1848.0, 1848.0, 1605.0,  596.0, &
     596.0,  333.0,  333.0,  333.0,  899.0, &
    5660.0, 5660.0, 2649.0, 2298.0, 2298.0, &
    2298.0,  209.0,  982.0,  309.0,  309.0, &
     309.0,  309.0,  589.0,  589.0, 5969.0, &
     589.0,  589.0,  589.0,  677.0, 1422.0, &
     589.0,  589.0,  589.0, 1422.0, 1422.0, &
    1422.0, 1028.0,  589.0/
!-----------
!  CHECK FOR DEBUG.
!-----------
CALL DBCHK (DEBUG,'DGF',3,ICYC)
IF(DEBUG) WRITE(JOSTND,5)ICYC
5 FORMAT(' ENTERING SUBROUTINE DGF  CYCLE =',I5)
!-----------
!  INITIALIZE VARIABLES FOR DBH >= 5 IN
!-----------
SDQGE5=0.
T=0.
BAGE5=0.
QMDGE5=0.
!----------
!  COMPUTE BASAL AREA IN TREES >= 5.0" DBH,
!  AND QMD FOR TREES >= 5.0" DBH
!----------
DO 20 ISPC=1,MAXSP
I1=ISCT(ISPC,1)
IF(I1 .LE. 0)GO TO 20
I2=ISCT(ISPC,2)
DO 10 I3=I1,I2
I=IND1(I3)
P=PROB(I)
D=DIAM(I)
IF(D.LT.5.)GO TO 10
SDQGE5=SDQGE5+P*(D)**2.
T=T+P
BAGE5=BAGE5+D*D*P*0.005454154
10 CONTINUE
20 CONTINUE
IF(SDQGE5.GT.0.)QMDGE5=SQRT(SDQGE5/T)
IF(BAGE5 .LE. 0.) BAGE5 = 10.
!----------
!  BEGIN SPECIES LOOP.
!----------
DO 200 ISPC=1,MAXSP
I1=ISCT(ISPC,1)
IF(I1.EQ.0) GO TO 200
I2=ISCT(ISPC,2)
CONSPP= DGCON(ISPC) + COR(ISPC)
!----------
!  BEGIN TREE LOOP WITHIN SPECIES ISPC.
!----------
DO 100 I3=I1,I2
I=IND1(I3)
D=DIAM(I)
CR=ICR(I)
IF(CR .LE. 0.) CR= 10.
BAL = (1.0 - (PCT(I)/100.)) * BA
!
IF (D.LE.0.0) GO TO 100
!----------
!  BOUND QMD FOR SOME SPECIES
!----------
SELECT CASE (ISPC)
  CASE(1,2,10,37:39,59)                        ! SC,JP,TA,BH,PH,SH,BG
    IF(QMDGE5 .GT. 13.) QMDGE5=13.
  CASE(11)                                     ! WC
    IF(QMDGE5 .GT. 15.) QMDGE5=15.
  CASE(17,30:33)                               ! EC,WO,SW,BR,CK
    IF(QMDGE5 .GT. 25.) QMDGE5=25.
  CASE DEFAULT
END SELECT
!----------
!  BOUND CROWN RATIO FOR SOME SPECIES
!----------
SELECT CASE (ISPC)
  CASE(17)                                     !EC
    IF(CR .GT. 60.) CR=60.
  CASE(60,65)                                  !SY,BL
    IF(CR .GT. 85.) CR=85.
  CASE DEFAULT
END SELECT
!----------
!  RELATIVE DBH
!----------
RELDBH = 0.0
RELDBHSQ = 0.0
IF(QMDGE5.GT.0.0) RELDBH=D/QMDGE5
IF(QMDGE5.GT.0.0) RELDBHSQ=D*D/QMDGE5
!
IF(DEBUG) WRITE(JOSTND,*) 'IN DGF-I= ',I,' D=',D, &
   ' RELDBH= ',RELDBH,' SI= ',SITEAR(ISPC), &
   'BAGE5=',BAGE5,'CR=',CR,' QMDGE5= ',QMDGE5, &
   ' BAL= ',BAL,' CONSPP= ',CONSPP
!----------
!  CALCULATION OF DDS FOR LAKE STATES VARIANT
!----------
DDS= CONSPP &
      + INTERC(ISPC) &
      + VDBHC(ISPC) * 1./D &
      + DBHC(ISPC) * D &
      + DBH2C(ISPC) * D * D &
      + RDBHC(ISPC) * RELDBH &
      + RDBHSQC(ISPC) * RELDBHSQ &
      + CRWNC(ISPC) * CR &
      + CRSQC(ISPC) * CR*CR &
      + SBAC(ISPC) * BAGE5 &
      + BALC(ISPC) * BAL &
      + SITEC(ISPC) * SITEAR(ISPC)
!----------
!  IF DEBUGGING, OUTPUT TERMS IN THE EQUATION
!----------
IF(DEBUG) THEN
  D1= INTERC(ISPC)
  D2= VDBHC(ISPC) * 1./D
  D3= DBHC(ISPC) * D
  D4= DBH2C(ISPC) * D * D
  D5= RDBHC(ISPC) * RELDBH
  D6= RDBHSQC(ISPC) * RELDBHSQ
  D7= CRWNC(ISPC) * CR
  D8= CRSQC(ISPC) * CR*CR
  D9= SBAC(ISPC) * BAGE5
  D10=BALC(ISPC) * BAL
  D11=SITEC(ISPC) * SITEAR(ISPC)
!
  WRITE(JOSTND,9090)D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,DDS
9090   FORMAT(' IN DGF ,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,DDSO=', &
     /,15F9.3)
ENDIF
!----------
!  CALCULATION OF DDS FOR LAKE STATES VARIANT WAS BUILT
!  DDS OUTSIDE BARK, CONVERT TO INSIDE BARK
!----------
IF(DDS.LT.-9.21) DDS=-9.21
DIAGRO= SQRT(DBH(I)*DBH(I)+EXP(DDS))-DBH(I)
BARK=BRATIO(ISPC,DBH(I),HT(I))
DIAGRI= DIAGRO*BARK
DDS=ALOG(((DBH(I)*BARK+DIAGRI)**2.0)-(DBH(I)*BARK)**2.0)
IF(DEBUG) WRITE(JOSTND,*) 'IN DGF-I= ',I,' D=',D, &
   ' DIAGRO= ',DIAGRO,' DIAGRI= ',DIAGRI,' BARK= ',BARK, &
   ' DDS= ',DDS
!
IF(DDS.LT.-9.21) DDS=-9.21
WK2(I)=DDS
!----------
!  END OF TREE LOOP.  PRINT DEBUG INFO IF DESIRED.
!----------
IF(DEBUG)THEN
  WRITE(JOSTND,9095) I,ISPC,DDS
9095   FORMAT(' IN DGF, I=',I4,',  ISPC=',I3,',  DDS=',F7.4)
  WRITE(JOSTND,*)' IN DGF, DG= ',(DG(IK), IK= 1,50)
ENDIF
!
100 CONTINUE
!----------
!  END OF SPECIES LOOP.
!----------
200 CONTINUE
IF(DEBUG) WRITE(JOSTND,9100)ICYC
9100 FORMAT(' LEAVING SUBROUTINE DGF  CYCLE =',I5)
!
RETURN
!
!
ENTRY DGCONS
!----------
!  ENTRY POINT FOR LOADING COEFFICIENTS OF THE DIAMETER INCREMENT
!  MODEL THAT ARE SITE SPECIFIC AND NEED ONLY BE RESOLVED ONCE.
!  OBSERV CONTAINS THE NUMBER OF
!  OBSERVATIONS BY SPECIES FOR THE UNDERLYING MODEL
!  (THIS DATA IS USED BY **DGDRIV** FOR CALIBRATION).
!
!----------
!  ENTER LOOP TO LOAD SPECIES DEPENDENT VECTORS.
!----------
DO ISPC=1,MAXSP
  DGCON(ISPC)=0.
  ATTEN(ISPC)=OBSERV(ISPC)
  SMCON(ISPC)=0.
!----------
!  IF READCORD OR REUSCORD WAS SPECIFIED (LDCOR2 IS TRUE) ADD
!  LN(COR2) TO THE BAI MODEL CONSTANT TERM (DGCON).  COR2 IS
!  INITIALIZED TO 1.0 IN BLKDATA.
!----------
  IF (LDCOR2.AND.COR2(ISPC).GT.0.0) DGCON(ISPC)=DGCON(ISPC) &
      + ALOG(COR2(ISPC))
ENDDO
!
RETURN
END
