!== last modified  10-4-2002
  SUBROUTINE SF_HS(JSP,GEOSUB,IEXTRA,SETOPT,DBH,TOTALH,DBTBH, &
                   HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,DIB,H)

!     FIND THE HEIGHT AT WHICH A SPECIFIED DIB OCCURS.
!      (0 IS RETURNED IF THE DIB IS IMPLAUSIBLE)
!     CONVERGENCE REQUIRES THESE TWO CONDITIONS:
!      1. TOLERANCE ON THE RETURNED HEIGHT IS PROBABLY BETTER
!         THAN (TOL * TOTAL HEIGHT).
!      2. PREDICTED DIB IS WITHIN EPSILON OF REQUESTED DIB
!         ON THE NEXT TO LAST ITERATION.
!
!     DIB   input    R*4   Diameter inside bark
!     H     output   R*4   Section height
!                             Warnings are written on unit 6 if this value is suspect.
!                             IF H=0 on return   NO SOLUTION was found
!
  CHARACTER*2 GEOSUB
  INTEGER INEEDSL,JSP,IEXTRA,SETOPT(6),IBREAK
  INTEGER ITER
  REAL hextra(2), z(2), slope
  REAL RHI1, RHI2, RHC, RHLONGI, F
  REAL RHFW(4),RFLW(6),TAPCOE(12),FMOD(3),PINV_Z(2)
  REAL DBH, TOTALH, DBTBH,RH,TOOHIGH,TOOLOW,DBT,DOB
!
  REAL*4 DIB,H,HI2,DI2,HI1,DI1,DBASE,RZ,D,ERR,ADJUST,EHIGH
  REAL*4 HHIGH,HLOW,HTRY,DTRY,D1,D2,ELOW,EPS,EPSILON,TOL

  DATA EPSILON /0.001/
  DATA TOL /0.0005/

!              Step 1. Calculate diameter at inflection point (DI2)
!                  and inflection height (HI2)
  RHI1 = RHFW(1)
  RHI2 = RHFW(2)
  RHC = RHFW(3)
  RHLONGI = RHFW(4)

  INEEDSL = 0
  RH = H/TOTALH
  HI2 = RHI2*TOTALH
  TOOHIGH = TOTALH
  TOOLOW = 0.
  CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH,HI2, &
             HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,DI2)
  IF(JSP.GE.22 .AND. JSP.LE.30) THEN
    CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HI2,DI2,DOB,DBT)
  ENDIF
!
!              Make a first guess for H; choice of formula
!              depends on whether we are above or below inflection.
  IF(DIB.GT.DI2) THEN
!                                          check the straight region
    TOOHIGH = HI2
    if(RHLONGI.gt.0.) then
      HI1=RHI1*TOTALH
      CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH, &
                 HI1,HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,DI1)
      IF(JSP.GE.22 .AND. JSP.LE.30) THEN
         CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HI1,DI1,DOB,DBT)
      ENDIF
      if (DIB.lt. DI1) then
            TOOLOW=DI1
            RH = RHI2 - (RHI2-RHI1)*(DIB-DI2)/(DI1-DI2)
            GO TO 100
      ENDIF
      TOOHIGH=DI1
    else
      DI1=DI2
    ENDIF
!                                 SOLUTION IS BELOW THE INFLECTION

    CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH, &
               0.0,HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,DBASE)
    IF(JSP.GE.22 .AND. JSP.LE.30) THEN
       CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,0.0,DBASE,DOB,DBT)
    ENDIF
    IF(DBASE.LE.DIB) THEN
      H=0.0
      INEEDSL=0
      RETURN
    ENDIF
!                Note: if there were an analytic solution, it would go here
    RZ= ((DIB-DI1)/(DBASE-DI1))**0.25
    RH= (1.0-RZ)*RHI1
  ELSE
!                                 SOLUTION IS ABOVE THE INFLECTION
    TOOLOW=DI2
    RZ =1.0-(DIB/DI2)**2.0
    RH = RHI2 + (1-RHI2)*RZ
  ENDIF
!
100 H=RH*TOTALH
  INEEDSL=1
  IBREAK=0
!
!  NEWTON SOLUTION FOR ROOTS
!
110 ITER=0
20 ITER=ITER+1
  IF(ITER.GT.30) THEN
!          WRITE(6,21) DIB,TOTALH,H-ADJUST,ERR
!   21     FORMAT(/' INITIAL FAILURE TO CONVERGE IN ROUTINE SF_HS',
!     &           /'  DIB, TOTALH, H, ERROR =',4F9.3  /
!     2           / ' WILL TRY A BISECTION APPROACH.')
      GO TO 200
  ENDIF
!
  CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH,H, &
             HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,D)
  IF(JSP.GE.22 .AND. JSP.LE.30) THEN
      CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HI2,D,DOB,DBT)
  ENDIF
  ERR= D - DIB
  IF (ERR.LT. 0.0) TOOHIGH=H
  ADJUST= - ERR/SLOPE
  H=H+ADJUST
  IF(H.GT.TOTALH) H=(H-ADJUST+TOTALH)/2.
  IF(H.LT.0.0) H=(H-ADJUST)/2.
  RH=H/TOTALH
  IF(ABS(ADJUST).GT.TOL*TOTALH .OR. ABS(ERR).GT.EPSILON) GO TO 20

  IF (SLOPE.gt. 0.) then
      IF(IBREAK.lt.2) then
           IBREAK=IBREAK+1
           if(IBREAK.eq.1)  H=0.8*H
           if(IBREAK.eq.2) H = H +(TOTALH-H)*0.25
           GO TO 110
         ENDIF
!          WRITE(6,22) DBH, TOTALH, DIB, H
!22        format(' SF_HS is probably returning the wrong root.'/
!     1           '     Program must be enhanced to guard against this.'/
!     2           ' DBH, TOTAL HT, (D, H) =',4f9.3 )
    ENDIF

!
!  NOTE: 1 STEP HAS BEEN TAKEN AFTER CONVERGENCE CRITERIA MET.
!
  INEEDSL=0
  RETURN
!                     MODIFIED   BISECTION APPROACH
!                       (slightly easier convergence criteria)
200 INEEDSL=0
  ITER=0
  HHIGH=TOOHIGH
  HLOW=TOOLOW
  CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH, &
             HHIGH,HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,D1)
  IF(JSP.GE.22 .AND. JSP.LE.30) THEN
      CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HHIGH,D1,DOB,DBT)
  ENDIF
  CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH,HLOW, &
             HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,D2)
  IF(JSP.GE.22 .AND. JSP.LE.30) THEN
      CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HLOW,D2,DOB,DBT)
  ENDIF
  EHIGH=D1-DIB
  ELOW=D2-DIB
  if(ELOW * EHIGH .gt. 0.) then
!               write(6,201)
!201            FORMAT(' SF_HS  BISECTION APPROACH FAILED.'/
!     1' Either a programming error or a major problem with user input.')
           H=0
           RETURN
        ENDIF


  EPS=EPSILON*2.0
210 HTRY = (HHIGH+HLOW)/2.

  CALL SF_DS(JSP,geosub,IEXTRA,SETOPT,ineedsl,slope,DBH,TOTALH,HTRY, &
             HEXTRA,Z,RHFW,RFLW,TAPCOE,F,FMOD,PINV_Z,DTRY)
  IF(JSP.GE.22 .AND. JSP.LE.30) THEN
      CALL BRK_UP(JSP,GEOSUB,DBH,TOTALH,DBTBH,HTRY,DTRY,DOB,DBT)
  ENDIF
  ERR = DTRY-DIB
  IF  (ABS(ERR).LT.EPS) then
          H=HTRY
          RETURN
         ENDIF
  if(ITER.gt. 40) then
!          write(6,211)
!211       format(' SF_HS   Bisection approach failed to converge.')
      H=HTRY
      RETURN
    ENDIF
  ITER=ITER+1
  if(ERR.gt. 0.0) then
        HLOW=HTRY
      else
        HHIGH=HTRY
      ENDIF
   GO TO 210

  END

