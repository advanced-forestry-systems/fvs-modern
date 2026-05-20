!== last modified  05-16-2013
! YW 2022/0808 added STUMP to call DVEST subroutine
  SUBROUTINE GROSSVOL(REGN,FORST,VOLEQ,MTOPP,MTOPS,STUMP,DBHOB, &
       DRCOB,HTTYPE,HTTOT,HTLOG,HT1PRD,HT2PRD,UPSHT1,UPSHT2,UPSD1, &
       UPSD2,HTREF,AVGZ1,AVGZ2,FCLASS,DBTBH,BTR,VOL,LOGVOL,LOGDIA, &
       LOGLEN,TLOGS,NOLOGP,NOLOGS,CUTFLG,BFPFLG,CUPFLG,CDPFLG, &
       SPFLG,CONSPEC,PROD,HTTFLL,live,BA,SI,CTYPE,ERRFLAG)
!
!--   THERE ARE CALLS TO THE FOLLOWING SUBROUTINES
!--     PROFILE - EXTERNAL
!--     DVE     - EXTERNAL
!--     R4VOL  - EXTERNAL
!--     R6VOL  - EXTERNAL
!--     BLMVOL - EXTERNAL
!--     R8VOL  - EXTERNAL
!--     R10VOL - EXTERNAL
!--     R12VOL - EXTERNAL
!******DECLARE VARIABLES

  CHARACTER*1 HTTYPE,LIVE,CTYPE
  CHARACTER*2 FORST,PROD
  CHARACTER*3 MDL
  character*4 conspec
  CHARACTER*10 VOLEQ

!     MERCH VARIABLES
  INTEGER REGN,HTTFLL,BA,SI
  REAL STUMP,MTOPP,MTOPS
  INTEGER CUTFLG,BFPFLG,CUPFLG,CDPFLG,SPFLG,ERRFLAG
!     TREE VARIABLES
  REAL HTTOT,HT1PRD,HT2PRD,THT1
  REAL DBHOB,DRCOB,DBTBH,BTR
  INTEGER FCLASS,HTLOG
!     3RD POINT VARIABLES
  REAL UPSD1,UPSD2,UPSHT1,UPSHT2,AVGZ1,AVGZ2
  INTEGER HTREF,I,J
!     OUTPUTS
  REAL VOL(15),LOGVOL(7,20),NOLOGP,NOLOGS
  REAL LOGDIA(21,3),LOGLEN(20),BOLHT(21)
  INTEGER TLOGS

!     TLOGS =  0
  MDL = VOLEQ(4:6)
! When total height is entered, the height type has to be feet. (2013/05/16)
  IF(HTTOT.GT.0) HTTYPE = 'F'

  IF(MDL.EQ.'FW2' .OR. MDL.EQ.'fw2' .OR. MDL.EQ.'FW3' .OR. &
     MDL.EQ.'fw3' .OR. MDL.EQ.'CZ2' .OR. MDL.EQ.'cz2' .OR. &
     MDL.EQ.'CZ3' .OR. MDL.EQ.'cz3' .OR. MDL.EQ.'WO2' .OR. &
     MDL.EQ.'wo2' .OR. MDL.EQ.'F32' .OR. MDL.EQ.'f32' .OR. &
     MDL.EQ.'F33' .OR. MDL.EQ.'f33' .OR. MDL.EQ.'JB2' .OR. &
     MDL.EQ.'jb2') THEN
!***********************
!      INGY MODELS     *
!      REGION 2 MODELS *
!      REGION 5 MODELS *
!***********************

    CALL PROFILE (REGN,FORST,VOLEQ,MTOPP,MTOPS,STUMP,DBHOB,HTTYPE, &
         HTTOT,HTLOG,HT1PRD,HT2PRD,UPSHT1,UPSHT2,UPSD1,UPSD2,AVGZ1, &
         AVGZ2,HTREF,DBTBH,BTR,LOGDIA,BOLHT,LOGLEN,LOGVOL,VOL,TLOGS, &
         NOLOGP,NOLOGS,CUTFLG,BFPFLG,CUPFLG,CDPFLG,SPFLG,DRCOB, &
         CTYPE,FCLASS,PROD,ERRFLAG)

  ELSEIF (MDL.EQ.'MAT' .OR. MDL.EQ.'mat') THEN
!***********************
!      REGION 4 MODEL  *
!***********************

    CALL R4VOL(REGN,VOLEQ,MTOPP,HTTOT,DBHOB,HT1PRD,VOL,NOLOGP, &
                 NOLOGS,LOGDIA,LOGLEN,LOGVOL,BOLHT, &
               CUTFLG,BFPFLG,CUPFLG,CDPFLG,SPFLG,ERRFLAG)
      TLOGS = ANINT(NOLOGP + NOLOGS)
  ELSEIF (MDL.EQ.'TRF' .OR. MDL.EQ.'trf')THEN
!********************************
!      PNW terif VOLUME EQUATION
!*******************************
    CALL PNWTARIF(VOLEQ,DBHOB,HTTOT,MTOPP,VOL,BFPFLG,CUPFLG, &
                      ERRFLAG)

  ELSEIF (VOLEQ(1:1).EQ.'6') THEN
!*******************************
!     REGION 6 VOLUME ROUTINES *
!*******************************
      tht1 = 0.0
      if(httot.gt.0) then
          tht1 = httot
      elseif(HTTYPE.EQ.'L' .OR. HTTYPE.EQ.'l') then
          tht1 = ht1prd
      else
          tht1 = ht1prd
      endif

      nologp=0.0
      nologs=0.0
!          do 66 i=1,20
!            lggrd(i)=' '
!   66     continue

! ***   If not top dib is specified then
! ***   use the one from the vol equation

      IF((MDL.EQ.'BEH' .OR. MDL.EQ.'beh').AND. &
                              (bfpflg.eq.1 .or. cupflg.eq.1)) then
         IF(MTOPP.EQ.0) MTOPP = 6.0
         CALL R6VOL(VOLEQ,FORST,DBHOB,BTR,FCLASS,MTOPP,THT1, &
              HTTYPE,VOL,LOGVOL,NOLOGP,LOGDIA,LOGLEN,DBTBH,HT1PRD, &
              CTYPE,ERRFLAG)
         TLOGS = ANINT(NOLOGP)
      ELSEIF((MDL.EQ.'DVE'.OR.MDL.EQ.'dve') .AND. BFPFLG.EQ.1) THEN
         CALL R6VOL2(VOLEQ,DBHOB,HTTOT,VOL,ERRFLAG)
      ELSE
        ERRFLAG = 1
      ENDIF


  ELSEIF(VOLEQ(1:1).EQ.'B' .or. voleq(1:1).eq.'b') THEN
!**************************
!     BLM VOLUME ROUTINES *
!**************************
      IF(fclass.eq.0) THEN
           ERRFLAG = 2
      ELSE
           call BLMVOL(VOLEQ,MTOPP,HTTOT,HT1PRD,DBHOB,HTTYPE,FCLASS, &
               VOL,LOGDIA,LOGLEN,LOGVOL,TLOGS,NOLOGP,NOLOGS, &
               BFPFLG,CUPFLG,ERRFLAG)

      ENDIF

  ELSEIF (VOLEQ(1:1).EQ.'8') THEN
!*********************
!      REGION 8 MODEL  *
!*********************

      CALL R8VOL (VOLEQ,DBHOB,HTTOT,UPSHT1,HT1PRD,MTOPP,PROD,VOL, &
                  FORST,SI,BA,CTYPE,BFPFLG,CUPFLG,SPFLG,ERRFLAG)

  ELSEIF (MDL.EQ.'DEM' .OR. MDL.EQ.'dem' .OR. &
          MDL.EQ.'CUR' .OR. MDL.EQ.'cur' .OR. &
          MDL.EQ.'BRU' .OR. MDL.EQ.'bru') THEN
!************************
!      REGION 10 MODEL  *
!************************
     IF (VOLEQ(1:3).EQ.'A01'.OR.VOLEQ(1:3).EQ.'A02' .or. &
          VOLEQ(1:3).EQ.'a01'.OR.VOLEQ(1:3).EQ.'a02' ) THEN

        CALL R10VOL(VOLEQ,MTOPP,MTOPS,HTTOT,HT1PRD,DBHOB, &
             HTTYPE,VOL,NOLOGP,NOLOGS,TLOGS,LOGLEN,LOGDIA,LOGVOL, &
             BFPFLG,CUPFLG,SPFLG,ERRFLAG)
     ELSE

        CALL PROFILE(REGN,FORST,VOLEQ,MTOPP,MTOPS,STUMP,DBHOB, &
           HTTYPE,HTTOT,HTLOG,HT1PRD,HT2PRD,UPSHT1,UPSHT2,UPSD1, &
           UPSD2,AVGZ1,AVGZ2,HTREF,DBTBH,BTR,LOGDIA,BOLHT,LOGLEN, &
           LOGVOL,VOL,TLOGS,NOLOGP,NOLOGS,CUTFLG,BFPFLG,CUPFLG, &
           CDPFLG,SPFLG,DRCOB,CTYPE,FCLASS,PROD,ERRFLAG)

     ENDIF

  ELSEIF(MDL.EQ.'SN2') THEN
!*********************
!      HAWAII MODEL  *
!*********************
    if(FCLASS.eq.0) then
         ERRFLAG = 2
    else
      call R12VOL(VOLEQ,MTOPP,HT1PRD,DBHOB,HTTOT,VOL,NOLOGP, &
                  NOLOGS,FCLASS,CUTFLG,BFPFLG,CUPFLG,ERRFLAG)
      TLOGS = ANINT(NOLOGP + NOLOGS)
    endif

  ELSEIF (MDL.EQ.'DVE' .OR. MDL.EQ.'dve') THEN
!*********************
!    DVE MODELS FOR  *
!      REGION 1      *
!      REGION 2      *
!      REGION 3      *
!      REGION 5      *
!      REGION 9      *
!      ARMY BASE     *
!*********************

     CALL DVEST (VOLEQ,DBHOB,DRCOB,HTTOT,MTOPP,FCLASS,HTLOG,HT1PRD, &
          HT2PRD,FORST,BTR,VOL,CUTFLG,BFPFLG,CUPFLG,CDPFLG,SPFLG, &
          PROD,HTTYPE,HTTFLL,NOLOGP,LIVE,BA,SI,CTYPE,ERRFLAG,MTOPS, &
          STUMP)

  ELSE
!          ERROR MESSAGE
     ERRFLAG = 1
     DO 5 I=1,15
        VOL(I)=0.0
5 CONTINUE
     DO 10 I=1,21
        BOLHT(I)=0.0
10 CONTINUE
     DO 15 I=1,7
        DO 16 J=1,20
           LOGVOL(I,J)=0.0
16 CONTINUE
15 CONTINUE
     DO 17,I=1,3
        DO 18,J=1,21
           LOGDIA(J,I)=0.0
18 CONTINUE
17 CONTINUE
     DO 19, I=1,20
        LOGLEN(I) = 0.0
19 CONTINUE
  ENDIF

4000 RETURN
  END
!********************************************************************
!********************************************************************
!--  VARIABLES OF INTEREST ARE:
!--******************************************************************
!--  Following variables have been be passed to the model *******
!--            optional variables are identified
!
!--  MERCHANDIZING VARIABLES
!***************************
!--  REGION - INTEGER - Region number used to set Regional Merchandizing Rules
!--  COR - CHARACTER - Flag to indicate Scribner table or Scribner
!--                 factor volumes. "Y" = table volumes, "N" = factor volumes
!--  EVOD - INTEGER - allow even or odd segment lengths
!--         segment options 11-14 allow odd lengths by definition
!--        1 = odd segment lengths allowed
!--        2 = only even segment lengths will be allowed
!--  MAXLEN - REAL - Maximum segment length
!--  MINLEN - REAL - Minimum segment length
!--  MERCHL - REAL - Minimum length of primary product a tree must have
!--                  must have be merchantable
!--  **TOP DIB TO USE**
!--  MTOPP - REAL - BdFt, CuFt and Cord Wood merch top for primary product
!--  MTOPS - REAL - CuFt and Cord Wood merch top for secondary product
!
!--  STUMP - REAL - height of stump in feet or fractions thereof.
!
!--  INTEGERS - Indicating which volumes are desired
!--    BFPFLG - Board foot primary product
!--    CUPFLG - Cubic foot primary product
!--    CDPFLG - Cord wood primary product
!--    SPFLG - Secondary product

!****** End of Merchandising parameters.
!
!    TREE CHARACTERISTICS
!************************
!--  VOLEQ - INTEGER - The volume equation number for this tree
!--  LIVE - CHAR -  Live/Dead code.  L=Live tree, D=Dead tree
!--  PROD - CHAR - Product code.  01 = sawtimber, 03 = pulp tree
!
!--  DBHOB - REAL - Diameter Breast Height Outside Bark
!--  DRCOB - REAL - Diameter Root Collar Outside Bark
!
!--  HTTYPE - CHAR - Tree height measured in feet (F) or logs (L)
!--  HTTOT - REAL - Tree height from the ground to the tip
!--  LOGHT - REAL - Tree height in number of logs to MTOPP
!--  HT1PRD - REAL - Tree height in feet to MTOPP
!--  HT2PRD - REAL - Tree height in feet to MTOPS
!
!--  Note measurements of bark and form are optional entries
!--  DBTBH - REAL - Double bark thickness at breast height
!--  FCLASS - INTEGER - Girard form class
!--  BTR - INTEGER - Bark thickness ratio
!--  UPSHT1 - REAL - The height of the first upper stem form measurement
!--  UPSHT2 - REAL - The height of the second upper stem form measurement
!--  UPSD1 - REAL - The measured diameter outisde bark at point UPSHT1.
!--  UPSD2 - REAL - The measured diameter outisde bark at point UPSHT2 .
!--  HTREF - REAL - The percent reference height for AVGZ1
!--  AVGZ1 - REAL - An average Z correction factor to use for point UPSHT1.
!--  AVGZ2 - REAL - An average Z correction factor to use for point UPSHT2.
!
!--  PCTDF - REAL - THE INTEGER VALUE REPRESENTING THE % OF DEFECT
!--          IN THE MAIN STEM PORTION OF A TREE - THIS DEFECT %
!--          WILL BE APPLIED TO BOTH BDFT AND CUFT PRODUCTS.
!--  PCTDF2 - REAL - THE INTEGER VALUE REPRESENTING THE % OF DEFECT
!--          IN THE TOPWOOD PORTION OF A TREE - THIS DEFECT %
!--          WILL BE APPLIED TO CUFT PRODUCTS.
!
!--  RETURNED VARIABLES
!--********************
!--  Following variables will be returned to the calling routine ******
!--             when requested
!--  VOL - REAL
!--    VOL(1) - Total volume for the tree in cubic feet, ground to tip.
!--              Does not include the limbs and roots.
!--    VOL(2) - Gross amount of board foot as primary product - Scribner.
!--    VOL(3) - Net amount of board foot as primary product - Scribner.
!--    VOL(4) - Gross amount of cubic foot as primary product - Smalian.
!--    VOL(5) - Net amount of cubic foot as primary product - Smalian.
!--    VOL(6) - Amount of cord wood as primary product.
!--    VOL(7) - Gross amount of cubic foot as secondary product - Smalian.
!--    VOL(8) - Net amount of cubic foot as secondary product - Smalian.
!--    VOL(9) - Amount of cord wood as secondary product.
!--    VOL(10) - Gross amount of board foot as primary product - INT'L 1/4"
!--    VOL(11) - Net amount of board foot as primary product - INT'L 1/4"
!--
!--  LOGLEN - REAL(20) - Log segment lengths computed by
!--                      subroutine "SEGMENT".  Does not include trim.
!--  LOGVOL(7,20) - REAL - The product volume of a log
!--    LOGVOL(1,X) - GROSS BOARD FT LOG VOLUME (20 LOGS)
!--    LOGVOL(2,X) - GROSS REMOVED BOARD FT LOG VOLUME (20 LOGS)
!--    LOGVOL(3,X) - NET BOARD FT LOG VOLUME (20 LOGS)
!--    LOGVOL(4,X) - GROSS CUBIC FT LOG VOLUME (20 LOGS)
!--    LOGVOL(5,X) - GROSS REMOVED CUBIC FT LOG VOLUME (20 LOGS)
!--    LOGVOL(6,X) - NET CUBIC FT LOG VOLUME (20 LOGS)
!--    LOGVOL(7,X) - GROSS INT'L 1/4 VOLUME (20 LOGS)
!
!--  LOGDIA - REAL(21,3) - Log end diameters
!       LOGDIA(x,1) = scaling dib
!       LOGDIA(x,2) = actual dib
!       LOGDIA(x,3) = actual dob
!
!--  BOLHT - REAL(21) - point on the bole where diameters were estimated
!--  *** Note LOGDIB and BOLHT are dimensioned 1 larger then the number of
!--     logs.  This is to hold the values for the small end of the top log.
!
!--  NOLOGP - integer - Average number of 16 ft logs, primary product.
!--  NOLOGS - integer - Average number of 16 ft logs, secondary product.
!
!--******************************************************************
!
!--  VARIABLES OF INTEREST USED WITHIN THIS ROUTINE ARE:
!
!--  LENMS - REAL - Length of primary product if a secondary product
!--                 is desired.
!--  LMERCH - REAL - Total merchantable length of stem returned
!--             from subroutine "MERLEN" and the merchantable
!--             length of stem with trim and unusable top
!--             removed returned from subroutine "SEGMNT"
!--  NUMSEG - INTEGER - The computed number of segments from
!--           subroutine "NUMSEG"
!--
!**************************************************************
!**************************************************************


