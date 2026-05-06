!== last modified  4-9-2002
  SUBROUTINE R2TAP(VOLEQ,DBHOB,D30,H30,TOP6,HT1,HT2,D2,DBTBH,errflg)
!       (TAPEQU-INT  DBHOB-R  HT1-R  HT2-R  D2-R)
!
!###########################################################
!      CHARACTER*3 TNIRP/'YES'/
!      CHARACTER*1 NEWPG
!###########################################################

  CHARACTER*3 MDL
  CHARACTER*10 VOLEQ

  REAL DBHOB, HT1, HT2, D2, D30, H30, TOP6, D30IB,DBTBH
  REAL*8 R2MBB(4,8), R2MBA(2,8), R2CZC(3,8)
  REAL*8 EDBH(2,8), R23PA(2,8), R23PB(2,8)
  REAL*8 DBH, HTTOT, HTUP, DIBUP, DIBCOR
  REAL*8 COFI, TERM1, TERM2, TERM3, TERM4
  REAL*8 E, R1, R2, R3, B1, B2, EDBHIB
  PARAMETER(E=2.7182818284)
  INTEGER SP,errflg

!***********************************************************
!*****   COEFFICIENTS FOR TWO PARAMETER MODEL   ************
!***********************************************************
!
!       COFFICIENTS FOR MAX AND BURKHART MODEL (B1,B2,B3,B4)
!  ASPEN (1)
  DATA ((R2MBB(I,J),I=1,4),J=1,8) &
                 /-5.18995,2.57262,-3.85160,117.934, &
!  LODGEPOLE (2)
                  -3.65010,1.45492,-2.20082, 52.058, &
!  PONDEROSA PINE BLACK HILLS (3)
                  -2.59737,0.96927,-1.43195, 50.867, &
!  PONDEROSA PINE OTHER (4)
                  -3.80739,1.75784,-3.56366, 55.776, &
!  WHITE FIR (5)
                  -2.91187,1.26772,-3.76391, 58.596, &
!  SUBALPINE FIR (6)
                  -3.11638,1.46021,-2.63725,105.472, &
!  ENGELMANN SPRUCE (7)
                  -2.26300,0.92540,-0.80682,382.694, &
!  DOUGLAS FIR (8)
                  -5.86345,2.98778,-4.12919, 82.838/

!          COFFICIENTS FOR MAX AND BURKHART MODEL (A1,A2)
!  ASPEN  (1)
  DATA ((R2MBA(I,J),I=1,2),J=1,8) &
                            /  0.69, 0.09, &
!  LODGEPOLE (2)
                               0.77, 0.11, &
!  PONDEROSA PINE BLACK HILLS (3)
                               0.75, 0.11, &
!  PONDEROSA PINE OTHER (4)
                               0.62, 0.13, &
!  WHITE FIR (5)
                               0.50, 0.13, &
!  SUBALPINE FIR (6)
                               0.55, 0.09, &
!  ENGELMANN SPRUCE (7)
                               0.65, 0.05, &
!  DOUGLAS FIR (8)
                               0.72, 0.12/

!     COFFICIENTS FOR CZAPLEWSKI SECOND-STAGE
!          MODEL (C1,C2,C3)
  DATA ((R2CZC(I,J),I=1,3),J=1,8) / &
!  ASPEN  (1)
   1.0,0.0,0.0, &
!  LODGEPOLE (2)
   1.0876,-0.0080764,0.0, &
!  PONDEROSA PINE BLACK HILLS (3)
   1.1331,-0.0095335,0.0, &
!  PONDEROSA PINE OTHER (4)
   1.1251,-0.0082315,0.0, &
!  WHITE FIR (5)
   1.0,0.0,0.0, &
!  SUBALPINE FIR (6)
   1.0,0.0,0.0, &
!  ENGELMANN SPRUCE (7)
   1.1263,-0.0080396,0.0, &
!  DOUGLAS FIR (8)
   1.0,0.0,0.0/

!***********************************************
!*****   COFFICIENTS FOR THREE PARAMETER MODEL
!***********************************************

!     DBH INSIDE BARK COEFFICIENTS (A,B)
!  ASPEN  (1)
  DATA ((EDBH(I,J),I=1,2),J=1,8) &
                            /  -0.018759, 0.931816, &
!  LODGEPOLE (2)
                               -0.161248, 0.977231, &
!  PONDEROSA PINE BLACK HILLS (3)
                               -0.75702,  0.961789, &
!  PONDEROSA PINE OTHER (4)
                               -0.653970, 0.950056, &
!  WHITE FIR (5)
                               -1.078182, 0.969977, &
!  SUBALPINE FIR (6)
                               -0.330462, 0.961494, &
!  ENGELMANN SPRUCE (7)
                               -0.664774, 0.995785, &
!  DOUGLAS FIR (8)
                               -0.710803, 0.942278/

!     MAX AND BURKHART COEFFICIENTS (A1,A2)
!  ASPEN  (1)
  DATA ((R23PA(I,J),I=1,2),J=1,8) &
                            /  0.69, 0.09, &
!  LODGEPOLE (2)
                               0.77, 0.11, &
!  PONDEROSA PINE BLACK HILLS (3)
                               0.75, 0.11, &
!  PONDEROSA PINE OTHER (4)
                               0.62, 0.13, &
!  WHITE FIR (5)
                               0.50, 0.13, &
!  SUBALPINE FIR (6)
                               0.55, 0.09, &
!  ENGELMANN SPRUCE (7)
                               0.65, 0.05, &
!  DOUGLAS FIR (8)
                               0.72, 0.12/

!     MAX AND BURKHART COEFFICIENTS (B3,B4)
!  ASPEN  (1)
  DATA ((R23PB(I,J),I=1,2),J=1,8) &
                            /  -4.29592, 117.934, &
!  LODGEPOLE (2)
                               -3.00413,  52.058, &
!  PONDEROSA PINE BLACK HILLS (3)
                               -2.40572,  50.867, &
!  PONDEROSA PINE OTHER (4)
                               -4.45250,  55.776, &
!  WHITE FIR (5)
                               -4.42031,  58.596, &
!  SUBALPINE FIR (6)
                               -3.12030, 105.472, &
!  ENGELMANN SPRUCE (7)
                               -1.11960, 382.694, &
!  DOUGLAS FIR (8)
                               -6.62047,  82.838/

!*************************************************************
!     THIS IS THE START OF THE MAIN LOGIG
! ************************************************************

!     CONVERT VALUES TO DOUBLE PRECISION FOR CACULATIONS

  DBH = DBHOB
  HTTOT = HT1
  HTUP = HT2
  MDL = VOLEQ(4:6)
  IF(VOLEQ(8:10).EQ.'746')THEN
      SP = 1
  ELSEIF(VOLEQ(8:10).EQ.'108')THEN
      SP = 2
  ELSEIF(VOLEQ(8:10).eq.'122' .AND. VOLEQ(1:3).EQ.'203')THEN
      SP = 3
  ELSEIF(VOLEQ(8:10).eq.'122' .AND. VOLEQ(1:3).EQ.'200')THEN
      SP = 4
  ELSEIF(VOLEQ(8:10).eq.'015')THEN
      SP = 5
  ELSEIF(VOLEQ(8:10).eq.'019')THEN
      SP = 6
  ELSEIF(VOLEQ(8:10).eq.'093')THEN
      SP = 7
  ELSEIF(VOLEQ(8:10).eq.'202')THEN
      SP = 8
  else
     errflg = 1
     d2 = 0
     return
  ENDIF
!*************************************************************
!                      2-PARAMETER MODEL               *******
!*************************************************************

  IF(MDL.EQ.'CZ2' .OR. MDL.EQ.'cz2') THEN

!   MAX AND BURKHART PORTION TO CALCULATE DIB ******************

    IF (HTUP .GT. HTTOT) HTUP = HTTOT
    TERM1 = R2MBB(1,SP)*((HTUP/HTTOT)-1)
    TERM2 = R2MBB(2,SP)*(((HTUP**2)/(HTTOT**2))-1)
    IF ( (HTUP/HTTOT) .LT. R2MBA(1,SP)) THEN
      COFI = 1.0
    ELSE
      COFI = 0.0
    ENDIF

    TERM3 = ((R2MBA(1,SP)-(HTUP/HTTOT))**2)*COFI
    TERM3 = R2MBB(3,SP)*TERM3

    IF ( HTUP/HTTOT .LT. R2MBA(2,SP)) THEN
      COFI = 1.0
    ELSE
      COFI = 0.0
    ENDIF

    TERM4 = ((R2MBA(2,SP)-(HTUP/HTTOT))**2)*COFI
    TERM4 = R2MBB(4,SP)*TERM4

    IF (TERM1+TERM2+TERM3+TERM4 .LT. 0) THEN
      DIBUP = 0.0
    ELSE
      DIBUP = DBH*SQRT(TERM1+TERM2+TERM3+TERM4)
    ENDIF


!    CZAPLEWSKI PORTION ****************************************

!    CALCULATE DIB

    TERM1 = R2CZC(1,SP)
    TERM2 = R2CZC(2,SP)*DBH
    TERM3 = R2CZC(3,SP)*(HTUP**2)

    DIBCOR = DIBUP*(TERM1+TERM2+TERM3)
    D2 = DIBCOR

!****************************************************************
!                  USE NEW THREE PARAMETER MODEL     ************
!****************************************************************

  ELSEIF(MDL.EQ.'CZ3' .OR. MDL.EQ.'cz3') THEN

    IF (DBTBH.GT.0) THEN
      EDBHIB = DBH - DBTBH
    ELSE
      EDBHIB = EDBH(2,SP)*DBH + EDBH(1,SP)
    ENDIF

    D30IB = D30*(1.0-(1.0-EDBHIB/DBH)*(1.0/(2.0-D30/DBH)))

    R1 = HTUP/HTTOT
    R2 = 4.5/HTTOT
    R3 = H30/HTTOT

    B1 = ((EDBHIB/DBH)**2)*((R3**2)-1.) + ((D30IB/DBH)**2)* &
                                          (1.-R2**2)

    IF(R3.LT.R23PA(1,SP)) THEN
      B1 = B1 + R23PB(1,SP)*((R2**2)-1.)*((R23PA(1,SP)-R3)**2)
    ENDIF

    IF(R2.LT.R23PA(1,SP)) THEN
      B1 = B1 + R23PB(1,SP)*(1-(R3**2))*((R23PA(1,SP)-R2)**2)
    ENDIF

    IF(R3.LT.R23PA(2,SP)) THEN
      B1 = B1 + R23PB(2,SP)*((R2**2)-1.)*((R23PA(2,SP)-R3)**2)
    ENDIF

    IF(R2.LT.R23PA(2,SP)) THEN
      B1 = B1 + R23PB(2,SP)*(1.-(R3**2))*((R23PA(2,SP)-R2)**2)
    ENDIF
    B1 = B1/(((1.-R2)*(1.-R3**2))-((1.-R2**2)*(1.-R3)))
    B2 = B1*(R3-1.)-((D30IB/DBH)**2)

    IF(R3.LT.R23PA(1,SP)) THEN
      B2 = B2 + R23PB(1,SP)*((R23PA(1,SP)-R3)**2)
    ENDIF

    IF(R3.LT.R23PA(2,SP)) THEN
      B2 = B2 + R23PB(2,SP)*((R23PA(2,SP)-R3)**2)
    ENDIF
    B2 = B2/(1.-(R3**2))

    DIBCOR = B1 * (R1-1.)+B2*((R1**2)-1.)
    IF(R1.LT.R23PA(1,SP)) THEN
      DIBCOR = DIBCOR + R23PB(1,SP)*((R23PA(1,SP)-R1)**2)
    ENDIF

    IF (R1.LT.R23PA(2,SP)) THEN
      DIBCOR = DIBCOR + R23PB(2,SP)*((R23PA(2,SP)-R1)**2)
    ENDIF
    DIBCOR = DIBCOR *(DBH**2)

!********************************************
!  ABOVE A 6-INCH TOP LOGIC
!********************************************

    IF(TOP6.LE.0.0) THEN
      IF (DIBCOR.LE.0.0) THEN
         D2=0.0
      ELSE
         D2=SQRT(DIBCOR)
      ENDIF

    ELSE
      IF(HTUP.GE.TOP6) THEN
        IF (D30IB.LT.6.0) THEN
!              D2SQ = D30IB**2-(D30IB**2)*((HTUP-TOP6)/(HTTOT-TOP6))
!              D2 = SQRT(D2SQ)
          D2 = D30IB*((HTTOT-HTUP)/(HTTOT-TOP6))
        ELSE
!              D2SQ = 36.0 - (36.0*((HTUP-TOP6)/(HTTOT-TOP6)))
!              D2 = SQRT(D2SQ)
          D2 = 6.0*((HTTOT-HTUP)/(HTTOT-TOP6))
        ENDIF
      ELSE
        D2 = SQRT(DIBCOR)
      ENDIF
    ENDIF
  ENDIF

!*****************************************************************
!    THIS IS THE END OF THE MAIN LOGIC
!*****************************************************************


  RETURN
  END
