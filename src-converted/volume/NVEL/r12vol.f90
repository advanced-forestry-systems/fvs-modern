!== last modified  1-18-2013
! 01/18/2013 Added calculation for tip volume VOL(15)
! 2018/10/10 YW ADDED VALID SPECIES CHECK
! 2018/11/07 YW ADDED HTTOT TO R12VOL FOR TOTAL CUBIC VOLUME
  SUBROUTINE R12VOL(EQNUM,MTOPP,HT1PRD,DBHOB,HTTOT,VOL,NOLOGP, &
                    NOLOGS,FCLASS,CUTFLG,BFPFLG,CUPFLG,ERRFLAG)

!--   THIS SUBROUTINE DETERMINES THE VOLUME OF A TREE
!--   USING VOLUME DETERMINATION ROUTINES FOR HAWAII.
!--   THERE ARE CALLS TO THE FOLLOWING SUBROUTINES
!--       R12TAP
!**************************************************************

  character*10 eqnum

  REAL DBHOB, HT1PRD, LOGVOL(21),NOLOGP, NOLOGS, HT2, TCVOL
  REAL MTOPP, VOL(15), LVOL,HTTOT,TOPD

  INTEGER VTYPE,SEGNUM,FCLASS,CUTFLG,BFPFLG,CUPFLG,ERRFLAG,I,SP

!--   SET ALL POTENTIAL VOLUMES AND AVG NUM OF LOGS TO ZERO
  DO 10 I=1,15
      VOL(I) = 0.0
10 CONTINUE

  NOLOGP = 0.0
  NOLOGS = 0.0
  SEGNUM = 0
  ERRFLAG = 0

!--   IF DBHOB OR HT1PRD EQUALS ZERO THEN DON'T CALCULATE THE VOLUME

  IF ( DBHOB .LE. 1.0)THEN
     ERRFLAG = 3
     GOTO 1000
  ENDIF
  IF(HT1PRD .LE. 0.0 ) THEN
     ERRFLAG = 4
     GO TO 1000
  ENDIF
!--   CHECK VALID SPECIES (2018/10/10 YW)
  IF(EQNUM(8:10).EQ.'301')THEN
    SP = 1
  ELSEIF(EQNUM(8:10).EQ.'671')THEN
    SP = 2
  ELSEIF(EQNUM(8:10).EQ.'510' .AND. EQNUM(1:3).EQ.'H00')THEN
    SP = 3
  ELSEIF(EQNUM(8:10).EQ.'510' .AND. EQNUM(1:3).EQ.'H01')THEN
    SP = 4
  ELSE
    ERRFLAG = 6
    RETURN
  ENDIF

!*************************************************************
!--    DETERMINE TOTAL CUBIC FIBER CONTENT (ALL TREES)       *
!*************************************************************
!--   SUBROUTINE "R12TAP" IS THE PROFILE MODEL
  IF(CUTFLG.EQ.1) THEN
    HT2 = HT1PRD
    VTYPE = 1
!     TOTAL CUBIC VOLUME SHOULD USE HTTOT BASED ON SHARPNECK MODEL(YW 2018/11/07)
    TOPD = 0.0
    IF(HTTOT.LT.0.1)THEN
      HTTOT = HT1PRD
      TOPD = MTOPP
    ENDIF
!        CALL R12TAP(EQNUM,DBHOB,HT1PRD,HT2,MTOPP,FCLASS,VTYPE,TCVOL)
    CALL R12TAP(EQNUM,DBHOB,HTTOT,HT2,TOPD,FCLASS,VTYPE,TCVOL)
    VOL(1) = TCVOL
  ENDIF

!**************************************************************
!**************************************************************
!           BOARD FOOT CACULATIONS  (DBHOB GE 10.0)              *
!**************************************************************

!--   DETERMINE IF A BOARD FOOT PRODUCT REPRESTNTATION IS NEEDED

  IF(DBHOB .GE. 10.0) THEN
!--*************** USE STEM PROFILE MODEL TO DETERMINE VOLUME
      IF(BFPFLG.EQ.1) THEN
           VTYPE = 2
           SEGNUM = HT1PRD/8.15
           NOLOGP = SEGNUM
           DO 100, I=1,SEGNUM
                HT2 = I
                CALL R12TAP(EQNUM,DBHOB,HT1PRD,HT2,MTOPP,FCLASS, &
                                                  VTYPE,LVOL)
                VOL(10) = VOL(10) + LVOL
                LOGVOL(I) = LVOL
!                IF(I.EQ.2.OR.I.EQ.4) THEN
!                    TEMPVOL = LOGVOL(I) + LOGVOL(I-1)
!                ENDIF
100 CONTINUE
      ENDIF

!*****************
      IF (VOL(10) .LT. 0.0) THEN
           VOL(10) = 0.0
      ENDIF

!--   (ENDIF FOR BOARD FOOT EQUATIONS)
  ENDIF

!**************************************************************
!*      CUBIC FOOT MAIN STEM EQUATIONS  ALL TREES             *
!**************************************************************

  IF (CUPFLG.EQ.1) THEN
!--*************** USE STEM PROFILE MODEL TO DETERMINE VOLUME
      IF(DBHOB.GE.7.0) THEN
           VTYPE = 3
           SEGNUM = HT1PRD/8.15
           DO 200, I=1,SEGNUM
              HT2 = I
              CALL R12TAP(EQNUM,DBHOB,HT1PRD,HT2,MTOPP,FCLASS, &
                                                   VTYPE,LVOL)
              VOL(4) = VOL(4) + LVOL
200 CONTINUE
!--       (ENDIF FOR EQNUM(2) EQUAL TO A VALID EQUATION)
      ENDIF
!*************
      IF (VOL(4) .LT. 0.0) THEN
           VOL(4) = 0.0
      ENDIF
!--   (ENDIF FOR CUBIC FOOT EQUATIONS)
  ENDIF

!**************************************************************
!*      CORD WOOD MAIN STEM EQUATIONS  ALL TREES             *
!**************************************************************

!--   PUT LOGIC HERE FOR CORDS

!**************************************************************

    VOL(2) = VOL(10)
    VOL(3) = VOL(11)

!****************************************************************
!--  DETERMINE TOP WOOD CUBIC FOOT PRODUCTS (ALL TREES)
!--    ONLY ONE STEM PROFILE EQUATION MAY BE SPECIFIED FOR BDFT,
!--    CUFT, AND CORDS AND ONLY ONE TOP DIAMETER MAY BE SPECIFIED
!--    BEFORE TOP WOOD VOLUME CAN BE DETERMINED.
!****************************************************************


!***************************************************************

      IF (VOL(7) .LT. 0.0) THEN
         VOL(7) = 0.0
      ENDIF
!     Calculate stem tip volume
  IF(VOL(4).GT.0.0 .AND. VOL(1).GT.0.0) VOL(15)=VOL(1)-VOL(4)
  IF(VOL(15).LT.0.01) VOL(15)=0.0
1000 CONTINUE

  RETURN
  END

!**************************************************************
!**************************************************************

!--  VARIABLES OF INTEREST ARE:

!--  DBHOB - REAL - DIAMETER BREAST HEIGHT
!--  EQNUM - INTEGER - EQUATION NUMBER TO USE
!--    EQNUM(1) - BOARD FOOT MAIN STEM
!--    EQNUM(2) - CUBIC FOOT MAIN STEM
!--    EQNUM(3) - CORD WOOD MAIN STEM
!--    EQNUM(4) - CUBIC FOOT TOP WOOD
!--    EQNUM(5) - CORD WOOD TOP WOOD

!--  TOP DIB TO USE
!--    TOPDIA(1) - BOARD FOOT MAIN STEM
!--    TOPDIA(2) - CORD WOOD MAIN STEM
!--    TOPDIA(3) - CUBIC FOOT MAIN STEM
!--    TOPDIA(4) - CORD WOOD TOP WOOD
!--    TOPDIA(5) - CUBIC FOOT TOP WOOD

!********BOARD FOOT VOLUME EQUATION NUMBERS********
!--    THE FOLLOWING ARE REGION 5 STEM PROFILE MODELS
!      1251 = KOA
!      1252 = OHIA
!      1253 = ROBUSTA EUCALYPTUS
!      1254 = SALIGNA EUCALYPTUS

!********CUBIC FOOT VOLUME EQUATION NUMBERS********

!--    THE FOLLOWING ARE REGION 5 STEM PROFILE MODELS
!      1251 = KOA
!      1252 = OHIA
!      1253 = ROBUSTA EUCALYPTUS
!      1254 = SALIGNA EUCALYPTUS
!
!--  HT1PRD - REAL - TREE HEIGHT FROM GROUND TO TIP
!--  LENMS - REAL - LENGTH OF MAIN STEM IF A TOP WOOD PRODUCT
!--                 IS DESIRED
!--  MAXLEN - REAL - MAXIMUM SEGMENT LENGTH
!--  MINLEN - REAL - MINIMUM SEGMENT LENGTH
!--  MERCHL - INTEGER - MINIMUM PRODUCT LENGTH FOR A MERCH TREE
!--  NOLOGP - REAL - AVERAGE NUMBER OF 8 FOOT LOGS IN MAIN STEM
!--  NOLOGS - REAL - AVERAGE NUMBER OF 8 FOOT LOGS IN TOP WOOD
!--  SEGNUM - INTEGER - THE COMPUTED NUMBER OF SEGMENTS


!--  PCTDF - REAL - THE INTEGER VALUE REPRESENTING THE % OF DEFECT
!--          IN THE MAIN STEM PORTION OF A TREE - THIS DEFECT %
!--          WILL BE APPLIED TO BOTH BDFT AND CUFT PRODUCTS.
!--  STUMP - REAL - HEIGHT OF STUMP
!--  TAPEQU - REAL - TAPER EQUATION TO USE - FROM EQNUM(?)
!--  TRIM - REAL - TRIM LENGTH FOR EACH SEGMENT
!--  VOL - REAL
!--    VOL(1) -  TOTAL VOLUME FOR THE TREE IN CUBIC FOOT, DOES NOT
!--              INCLUDE LIMBS AND ROOTS.
!--    VOL(2) - GROSS AMOUNT OF BOARD FOOT PRODUCT IN MAIN STEM.
!--    VOL(3) - NET AMOUNT OF BOARD FOOT PRODUCT IN MAIN STEM.
!--    VOL(4) - GROSS AMOUNT OF CUBIC FOOT PRODUCT IN MAIN STEM.
!--    VOL(5) - NET AMOUNT OF CUBIC FOOT PRODUCT IN MAIN STEM.
!--    VOL(6) - AMOUNT OF CORD WOOD PRODUCT IN MAIN STEM
!--    VOL(7) - GROSS AMOUNT OF CUBIC PRODUCT IN TOP WOOD.
!--    VOL(8) - NET AMOUNT OF CUBIC PRODUCT IN TOP WOOD.
!--    VOL(9) - AMOUNT OF CORD WOOD PRODUCT IN TOP WOOD.
!--    VOL(10) - GROSS AMOUNT OF INT'L 1/4" BOARD FOOT PRODUCT IN MAIN STEM.
!--    VOL(11) - NET AMOUNT OF INT'L 1/4" BOARD FOOT PRODUCT IN MAIN STEM.
!--    VOL(12) - GROSS AMOUNT OF BOARD FOOT PRODUCT IN TOPWOOD.
!--    VOL(13) - NET AMOUNT OF BOARD FOOT PRODUCT IN TOPWOOD.
!--
!**************************************************************
