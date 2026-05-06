!== last modified  04-02-2014
  SUBROUTINE FWINIT(VOLEQ,DBHOB,HTTOT,MHT,MTOPP,UPSHT1,UPSHT2, &
       UPSD1,UPSD2,AVGZ1,AVGZ2,HTREF,DBTBH,JSP,RHFW,RFLW, &
       TAPCOE,F,SETOPT,NEXTRA,HEX,DEX,ZEX,PINV_Z,FMOD,btr,FCLASS, &
       ERRFLAG)
!**************************************************************
!        Initialiazes SF Subroutines using Flewelling profile models
!        Uses SF_TEST2 initializing code
!                  J. Flewelling  July, 1996
!                  Modified by K. Cormier Sept, 1996
!  04/02/2014 YW Changed UPSHT1 to 33.6 for 32 foot log when using FCLASS
!  07/12/2016 YW ADDED BTR VALUE FOR REGION 3 SANTA FE NF SPECIES 122,202,015,108
!  07/13/2023 YW Changed BTR for 301FW2W equatin for species 122,202,015,108 with R3 combined validation data
!     TREE VARIABLES
  REAL DBHOB,HTTOT,DBTBH,MHT,MTOPP,btr

!     3PT VARIABLES
  REAL UPSHT1,UPSHT2,UPSD1,UPSD2,AVGZ1,AVGZ2
  REAL*4 HEX(2),DEX(2),ZEX(2),FMOD(3),PINV_Z(2)
  REAL*4 RHFW(4),RFLW(6),TAPCOE(12),F

  INTEGER JSP, SETOPT(6),NEXTRA,htref,ERRFLAG,ICON,FCLASS,threeflag
!      LOGICAL METRIC
  Character*1 GEOCODE
!     character*1 unit
  Character*2 GEOSUB
  CHARACTER*3 SPEC
  CHARACTER*10 VOLEQ

!    Units in english or metric (E=english)
!      UNIT='E'
!                              use default settings
!    Apply constraints to 3-pt fit (1=yes/0=no)
  ICON=1
!    Stem form: user regin when available (1=yes/0=no)
!         IREGIONU(1)=1
!    Bark thickness: user regin when available (1=yes/0=no)
!         IREGIONU(2)=1
!
!    Upper stem measurement (1=DOB/2=DIB)
!      IDTYPE=1
!    Supplimental regression use intercept (1=yes/0=no)
!      IMOD_USE=1
!                Line from SF-SETERR
!      IER_UNIT=6
!                                                Other setup features
!                                                Note: Call changed from 1994.
!         FOLLOWING FROM SF_SET
!      if (unit.eq.'E' .or. unit.eq.'e') then
!          METRIC  = .false.
!          BH=4.5
!      else if  (unit.eq.'M' .or. unit.eq.'m') then
!          METRIC=.TRUE.
!          BH = 4.5 * 0.3048
!      else
!          IER=1
!      endif

!           set an upper bound on fractional changes;
!           predictions can exceed this only if measured data also exceeds.
  SETOPT(1)=ICON
!                              Region effects on bark predictions
  SETOPT(2)=1
  SETOPT(3)=1

!                          The following white noise option is ALWAYS set
!                            to zero (off). For future development.
  SETOPT(4)=0
!                           Numerical volume segments
!     SETOPT(5)= NVSEGU
!     SETOPT(6)= IVOLINT

!                       The following is the maximum allowed 3-pt fractional
!                         change in DBHOB. Unless a greater change is observed,
!                         this will become FMODMAX
!      SETC(1) = .15

!
!                               Prepare for modifier regressions
!      INT_USE = IMOD_USE
!     CALL SF_MODR0
!                                INGY
  ZEX(1) = 0
  ZEX(2) = 0
  HEX(1) = 0
  HEX(2) = 0
  DEX(1) = 0
  DEX(2) = 0
  NEXTRA = 0
  THREEFLAG = 0

!     MOVED THE FOLLOWING LINE TO END OF THE SETUP PHASE(YW 07/12/2016)
!      IF(BTR.GT.0.0 .AND. DBTBH.LE.0) DBTBH = DBHOB-(DBHOB*BTR/100.0)

  GEOCODE = VOLEQ(1:1)
  GEOSUB = VOLEQ(2:3)
  SPEC = VOLEQ(8:10)

!       ERROR CHECK THE GEOSUB CODE
  IF(GEOSUB.EQ.'OO' .OR. GEOSUB.EQ.'oo') GEOSUB = '00'

  IF(GEOCODE.EQ.'I' .OR. GEOCODE.EQ.'i' .OR. GEOCODE.EQ.'1') THEN
    IF (SPEC.EQ.'202'.or.spec.eq.'205'.or.spec.eq.'204') THEN
!     Douglas fir
      JSP = 11
    ELSEIF(SPEC.EQ.'073'.or.SPEC.EQ.'070') THEN
!     Western Larch
      JSP = 12
    ELSEIF(SPEC.EQ.'017') THEN
!     Grand fir
      JSP = 13
    ELSEIF(SPEC.EQ.'122') THEN
!     Ponderosa pine
      JSP = 14
    ELSEIF(SPEC.EQ.'108') THEN
!     Lodgepole pine
      JSP = 15
    ELSEIF(SPEC.EQ.'242'.or.SPEC.EQ.'240') THEN
!     Western Red Cedar
      JSP = 16
    ELSEIF(SPEC.EQ.'260'.or.SPEC.EQ.'263'.OR.SPEC.EQ.'264') THEN
!     Mountain Hemlock
      JSP = 17
    ELSEIF(SPEC.EQ.'119') THEN
!     White pine
      JSP = 18
    ELSEIF(SPEC.EQ.'093'.or.SPEC.EQ.'090') THEN
!     Engelmann Spruce
      JSP = 19
    ELSEIF(SPEC.EQ.'019') THEN
!     Subalpine fir
      JSP = 20
    ELSEIF(SPEC.EQ.'012') THEN
!     Balsam fir
      JSP = 21
    ELSE
      ERRFLAG = 1
      RETURN
    ENDIF

!                                WESTSIDE SPECIES
  ELSEIF(GEOCODE.EQ.'F' .OR. GEOCODE.EQ.'f') THEN
    IF (SPEC.EQ.'202'.or.spec.eq.'205'.or.spec.eq.'204') THEN
!     Douglas fir
      JSP = 3
    ELSEIF(SPEC.EQ.'263') THEN
!     Western Hemlock
      JSP = 4
    ELSEIF(SPEC.EQ.'242') THEN
!     Western Red Cedar
      JSP = 5
    ELSE
      ERRFLAG = 1
      RETURN
    ENDIF

!                                 REGION 2
  ELSE IF(GEOCODE.EQ.'2') THEN
    IF(SPEC.EQ.'122') THEN
!         Black Hills model
      IF(GEOSUB.EQ.'03') THEN
        JSP = 22
      ELSE
!         region wide and san juan
          JSP = 23
      ENDIF
    ELSEIF(SPEC.EQ.'108') THEN
!         Lodgepole model
       JSP=25
    ELSEIF(SPEC.EQ.'202') THEN
!         Douglas Fir model
       JSP = 26
    ELSEIF(SPEC.EQ.'015') THEN
!         White fir model
       JSP = 27
    ELSE IF(SPEC.EQ.'746')THEN
!         Aspen model
       JSP = 28
    ELSE
      ERRFLAG = 1
      RETURN
    ENDIF
!                                 REGION 4
  ELSE IF(GEOCODE.EQ.'4') THEN
    IF(GEOSUB.EQ.'07')THEN
      IF(SPEC.EQ.'093') THEN
!         Engelmann spruce model
        JSP = 24
      ELSE IF(SPEC.EQ.'122') THEN
!         R2 Ponderosa Pine model with R4 bark
        JSP = 23
      ELSE
        ERRFLAG = 1
        RETURN
      ENDIF
    ELSE
      ERRFLAG = 1
      RETURN
    ENDIF
!                                 REGION 3
!     REGION 3
  ELSE IF (GEOCODE.EQ.'3') THEN
     IF(GEOSUB.EQ.'00')THEN
       IF(SPEC.EQ.'122')THEN
         JSP = 29
       ELSEIF(SPEC.EQ.'202')THEN
          JSP = 26
          IF(BTR.LE.0) BTR = 88.85  !the BTR for combined R3 DF data
       ELSE
         ERRFLAG = 1
         RETURN
       ENDIF
!        GEOSUB 01 is for Santa Fe NF
     ELSEIF(GEOSUB.EQ.'01') THEN
        IF(SPEC.EQ.'122') THEN
          JSP = 29
          IF(BTR.LE.0) BTR = 89.12
!         SPECIES 202, 015, AND 108 USE REGION 2 PROFILE
        ELSEIF(SPEC.EQ.'108')THEN
          JSP = 25
          IF(BTR.LE.0) BTR = 93.26
        ELSEIF(SPEC.EQ.'202')THEN
          JSP = 26
          IF(BTR.LE.0) BTR = 89.72
        ELSEIF(SPEC.EQ.'015')THEN
          JSP = 27
          IF(BTR.LE.0) BTR = 91.16
        ELSE
          ERRFLAG = 1
          RETURN
        ENDIF
     ELSE
       ERRFLAG = 1
       RETURN
     ENDIF
!                                 REGION 10
  ELSE IF(GEOCODE.EQ.'A') THEN
      IF(SPEC.EQ.'042') THEN
!         Alaska yellow cedar
        JSP = 31
      ELSE IF(SPEC.EQ.'242') THEN
!         Western Red Cedar
        JSP = 32
      ELSE IF(SPEC.EQ.'098') THEN
!         Spruce
        IF(GEOSUB.EQ.'02') THEN
          JSP=35
        ELSE
          JSP = 33
        ENDIF
      ELSE IF(SPEC.EQ.'263'.OR.SPEC.EQ.'260'.OR.SPEC.EQ.'264') THEN
!         Hemlock
        IF(GEOSUB.EQ.'02') THEN
          JSP=36
        ELSE
          JSP = 34
        ENDIF
      ELSE
        ERRFLAG = 1
        RETURN
      ENDIF
  ELSE
    ERRFLAG = 1
    RETURN
  ENDIF

  IF(BTR.GT.0.0 .AND. DBTBH.LE.0) DBTBH = DBHOB-(DBHOB*BTR/100.0)
!-------------------------- end of set-up phase ------------------------------
!
!      WRITE(*,*)'FWINIT',JSP
!-----------------------------------------------------------------------------
!                 TWO POINT MODEL INITIALIZE,
!    check for total height
  IF(HTTOT.LE.0 .AND. MHT.GT.0) THEN
    CALL SF_2PTH(JSP,GEOSUB,SETOPT,DBHOB,MHT,MTOPP,DBTBH, &
           HTTOT,HEX,ZEX,F,RHFW,RFLW,TAPCOE)
    if(httot .le. 0) then   !added  10/4/2001
       errflag = 11
       return
    endif
    if(upsht1 .gt. 0) threeflag = 1
  ELSE
     CALL sf_2pt(JSP,GEOSUB,SETOPT,DBHOB,HTTOT,DBTBH,F,RHFW, &
                                                RFLW,TAPCOE)
  ENDIF
!          if (ierr.gt.0) then
!            write(6,106)  IERR
!106         format(' SF_2pt returned error code ', I6)
!            GOTO 1000
!          endif

!-----------------------------------------------------------------------------
!                                                    THREE POINT MODEL INITIALIZE
!                                                       Specify upper dob(s)
  IF(VOLEQ(6:6).EQ.'3'.OR.THREEFLAG.EQ.1)THEN
    IF (htref.gt.0)then
        NEXTRA = 1
        upsht1 = float(htref)*HTTOT/100.0
    ELSE IF(UPSHT1.GT.0.0) THEN
       IF(UPSHT2 .GT. 0.0) THEN
           NEXTRA = 2
       ELSE
           NEXTRA = 1
       ENDIF
    ELSE IF(FCLASS .GT. 0)THEN
       NEXTRA = 1
       UPSHT1 = 17.5
!          FCLASS for 32 foot log (YW 04/02/14)
       IF(VOLEQ(5:5).EQ.'3') UPSHT1 = 33.6
       UPSD1 = ANINT(DBHOB*FCLASS) / 100.0
    ENDIF

    IF (NEXTRA .GT. 0) THEN
      HEX(1) = UPSHT1
      HEX(2) = UPSHT2
!                                                      Specify an upper dib
      IF (AVGZ1 .ne. 0.0) THEN
        ZEX(1) = AVGZ1
        ZEX(2) = AVGZ2

        CALL SF_3z(JSP,GEOSUB,SETOPT,NEXTRA,DBHOB,HTTOT,HEX,ZEX, &
                         RFLW,RHFW,TAPCOE,F,PINV_Z,FMOD)

      ELSEIF(UPSD1 .GT. 0.0) THEN
        DEX(1) = UPSD1
        DEX(2) = UPSD2

        CALL sf_3pt(JSP,GEOSUB,setopt,DBHOB,HTTOT,DBTBH,NEXTRA, &
                    HEX,DEX,ZEX,RHFW,RFLW,TAPCOE,F,PINV_Z,FMOD)
!                                                       Specify upper z(s)
      ENDIF
    ELSE
      ERRFLAG = 9
    ENDIF
  ENDIF
!                    END OF INTIALIZATION ROUTINE
  RETURN
  END

! JSP VALUES
!
!    WESTSIDE:              3 - 5
!    INGY:                 11 - 21
!      Douglas fir            11
!      Western Larch          12
!      Grand fir              13
!      Ponderosa pine         14
!      Lodgepole pine         15
!      Western Red Cedar      16
!      Mountain Hemlock       17
!      White pine             18
!      Engelmann Spruce       19
!      Subalpine fir          20
!      Balsam fir             21
!
!    BLACK HILLS              22
!    SAN JUAN & R2 PP         23
!    DIXIE ES                 24
!    R2 LODGEPOLE             25
!    R2 DOUG FIR              26
!    R2 WHITE FIR             27
!    R10 ALASKA CEDAR         31
!    R10 WESTERN RED CEDAR    32
!
!
!

