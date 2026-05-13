!== last modified  09-15-2016 reconciled Vol. Eq. No. output from FVS with Cruise software-RNH
! 01/18/2013 added FIAVOLEQDEF, R5_PNWEQN and R6_PNWEQN for PNE FIA equations.
! 03/25/2014 changed default equation for Region 3 (R3_EQN) Ponderosa pine in the forest Apache Sitgreaves, Coconino, Kaibab and Tonto to 300FW2W122.
! 09/09/2016 Modified R3_EQN default equation for PP, DF, WF, and WP in Sante Fe NF
! 09/15/2016 Corrected R4 default equation for other species to DVE equation
! 04/12/2017 removed 532WO2W*** equation from R5_EQN per the email from Craig Bodenhausen (4/12/2017)
! 07/19/2021 Changed R8_CEQN to use the R8 new Clark equation 8*1CLKE***. Added 1 to TOPCODE array.
! 02/27/2025 Updated R6_EQN for species 81 Incense cedar to use I00FW2W073 in Willamette NF
! 2025/05/23 Updated R8_CEQN species list to include R8 cruise species
! 2024/09/09 Updated R3_EQN defaults for details see comment in R3_EQN
! 2025/11/24 Changed R6_EQN defaults for Red Alder in GP (3) and Siuslaw(12)
  SUBROUTINE VOLEQDEF (VAR,REGN,FORST,DIST,SPEC,PROD,VOLEQ,ERRFLAG)
!
!    SUBROUTINE WILL RETURN THE DEFAULT VOLUME EQUATION NUMBER
!        SPEC = 3 DIGIT FIA SPECIES CODE
  !DEC$ ATTRIBUTES DLLEXPORT::VOLEQDEF

  CHARACTER*2 FORST,PROD,DIST,VAR
  CHARACTER*10 VOLEQ
  INTEGER ERRFLAG,REGN,SPEC
!
!  The GETVARIANT routine is required by the NVEL system to find the
!  FVS variant. When VOLEQDEF is used by FVS, the call to VOLEQDEF always
!  carries the variant 2 character symbol so the GETVARIANT routine
!  is not used by FVS.
!
  IF(VAR .EQ. '  ') THEN
     CALL GETVARIANT(REGN,FORST,DIST,VAR)
  ENDIF
  IF(REGN .EQ. 1) THEN
     CALL R1_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.2)THEN
     CALL R2_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.3)THEN
     CALL R3_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.4)THEN
     CALL R4_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.5)THEN
     CALL R5_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.6)THEN
     CALL R6_EQN(VAR,FORST,DIST,SPEC,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.7)THEN
     CALL R7_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
!      ELSE IF((REGN.EQ.8).AND.(PROD.EQ.'02'))THEN
!         CALL R8_CEQN(FORST,DIST,SPEC,PROD,VAR,VOLEQ,ERRFLAG)
!      ELSE IF((REGN.EQ.8).AND.(PROD.EQ.'01'))THEN
!         CALL R8_BEQN(FORST,DIST,SPEC,PROD,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.8)THEN
     CALL R8_CEQN(FORST,DIST,SPEC,PROD,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.9)THEN
     CALL R9_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.10)THEN
     CALL R10_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  ENDIF

  RETURN
  END

!//////////////////////////////////////////////////////////////////
  SUBROUTINE GETVARIANT(REGN,FORST,DIST,VAR)
  CHARACTER*2 FORST,DIST,VAR
  CHARACTER VVER*7
  INTEGER REGN,FORNUM,DISTNUM

  IF(FORST(2:2) .LT. '0') THEN
    FORST(2:2) = FORST(1:1)
    FORST(1:1) = '0'
    IF(FORST(2:2) .LT. '0') FORST(2:2) = '0'
  ENDIF
  IF(DIST(2:2) .LT. '0') THEN
    DIST(2:2) = DIST(1:1)
    DIST(1:1) = '0'
    IF(DIST(2:2) .LT. '0') DIST(2:2) = '0'
  ENDIF
  READ(FORST,'(I2)')FORNUM
  READ(DIST,'(I2)')DISTNUM

  IF(REGN .EQ. 8)THEN
     VAR = 'SN'

  ELSE IF(REGN.EQ.1)THEN
     IF(FORNUM.EQ.4 .OR. FORNUM.EQ.5 .OR. FORNUM.EQ.17 .OR. &
        FORNUM.EQ.3 .OR. FORNUM.EQ.14 .OR. FORNUM.EQ.16)THEN
        VAR = 'IE'
     ELSE
        VAR = 'EM'
     ENDIF

  ELSE IF(REGN.EQ.5) THEN
!  DETERMINE VARIENT BY FOREST AND DISTRICT NUMBER
!  INLAND CALIFORNIA
     IF(FORNUM.EQ.5 .OR. FORNUM.EQ.6 .OR. FORNUM.EQ.8 .OR. &
        FORNUM.EQ.11 .OR. FORNUM.EQ.14) THEN
        VAR = 'CA'
!  SOUTHERN OREGON
     ELSEIF(FORNUM.EQ.9) THEN
        VAR = 'SO'
!  WESTERN SIERRA
     ELSEIF(FORNUM.EQ.17 .OR. FORNUM.EQ.16 .OR. FORNUM.EQ.15 .OR. &
            FORNUM.EQ.13 .OR. FORNUM.EQ.3) THEN
        VAR = 'WS'
!  KLAMATH/NORTHERN CALIFORNIA
     ELSEIF(FORNUM.EQ.5) THEN
        VAR = 'NC'
     ENDIF

  ELSE IF(REGN.EQ.6) THEN
!  DETERMINE VARIENT BY FOREST AND DISTRICT NUMBER
!  BLUE MTN VARIANT
     IF(FORNUM.EQ.4 .OR. FORNUM.EQ.7 .OR. FORNUM.EQ.14 .OR. &
        FORNUM.EQ.16) THEN
        VAR = 'BM'
!  EASTERN CASCADES
     ELSEIF(FORNUM.EQ.17 .OR. FORNUM.EQ.8 .OR. (FORNUM.EQ.3 .AND. &
            DISTNUM .EQ. 3) .OR. (FORNUM.EQ.6 .AND. &
           (DISTNUM.EQ.1.OR.DISTNUM.EQ.2.OR.DISTNUM.EQ.6))) THEN
!  MOUNT HOOD Barlow RD
        VAR = 'EC'
!  SOUTHERN OREGON
     ELSEIF(FORNUM.EQ.1 .OR. FORNUM.EQ. 2 .OR. FORNUM.EQ.20) THEN
        VAR = 'SO'
!  WESTERN CASCADES
     ELSEIF(FORNUM.EQ.5 .OR. FORNUM.EQ.15 .OR. FORNUM.EQ.18 .OR. &
            FORNUM.EQ.10 .OR. FORNUM.EQ.3 .OR. FORNUM.EQ.6) THEN
        VAR = 'WC'
!  PACFIC NORTHWEST
     ELSEIF(FORNUM.EQ.9 .OR. FORNUM.EQ.12) THEN
        VAR = 'PN'
!  NORTHERN CALIFORNIA
     ELSEIF(FORNUM.EQ.11) THEN
        VAR = 'NC'
     ELSEIF(FORNUM.EQ.21) THEN
        VAR = 'IE'
     ENDIF

  ELSE IF(REGN.EQ.7) THEN
!  DETERMINE VARIENT BY FOREST AND DISTRICT NUMBER
     IF(FORNUM.EQ.2) THEN
        VAR = 'WC'
     ELSEIF(FORNUM.EQ.3) THEN
        VAR = 'NC'
     ELSE
        VAR = 'SO'
     ENDIF

  ELSE IF (REGN .EQ. 9) THEN
     IF(FORNUM.EQ.13 .OR. FORNUM.EQ.10 .OR. FORNUM.EQ.3 .OR. &
          FORNUM.EQ.9 .OR. FORNUM.EQ.4 .OR. FORNUM.EQ.7 .OR. &
          FORNUM.EQ.2 .OR. FORNUM.EQ.6)THEN
        VAR = 'LS'
     ELSE IF(FORNUM.EQ.12 .OR. FORNUM.EQ.8 .OR. FORNUM.EQ.5)THEN
        VAR = 'CS'
     ELSE IF(FORNUM.EQ.21 .OR. FORNUM.EQ.20 .OR. FORNUM.EQ.19 .OR. &
             FORNUM.EQ.14 .OR. FORNUM.EQ.22)THEN
        VAR = 'NE'
     ENDIF

  ENDIF

  RETURN
  END
!
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R1_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
!      SEARCH BY FIA SPECIES CODE
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST,VAR
  INTEGER SPEC,ERRFLAG
  CHARACTER*10 EQNUM(81)
  INTEGER FIA(38), FIRST, HALF, LAST, DONE,FORNUM,I
! SPECIES
! Balsam fir,           White fir,           Grand fir,       Subalpine fir,   Western juniper,
! Rocky Mt Juniper,     Subalpine larch,     Western larch,   Engelmann spruce, White spruce,
! Blue spruce,          Whitebark pine,      Bristlecone pine,Pinyon Pine,      Lodgepole pine,
! Limber pine,          Western white pine,  Ponderosa pine,  Singleleaf pinyon,border pinyon,
! Douglas fir,          Pacific yew,         Western redcedar,Western hemlock,  Mountain hemlock,
! Other Softwoods,      Rocky Mountain maple,Paper birch,     curlleaf mtn-mahog,green ash,
! Cottonwood,           balsam poplar,       plains cottonwood,Quaking aspen,  black cottonwood,
! Narrowleaf cottonwood,Other hardwood,      Unknown

  DATA (FIA(I),I=1,38)/ 12,   15,   17,   19,   64, &
                        66,   72,   73,   93,   94, &
                        96,  101,  102,  106,  108, &
                       113,  119,  122,  133,  134, &
                       202,  231,  242,  263,  264, &
                       299,  321,  375,  475,  544, &
                       740,  741,  745,  746,  747, &
                       749,  998,  999/


  DATA (EQNUM(I),I=1,81)/ &
  'I00FW2W012','102DVEW017','I00FW2W017','I00FW2W019','102DVEW060', &
  '102DVEW106','I00FW2W019','I00FW2W073','I00FW2W093','102DVEW090', &
  '102DVEW090','I00FW2W119','102DVEW108','102DVEW106','I00FW2W108', &
  'I00FW2W073','I00FW2W119','I00FW2W122','102DVEW106','102DVEW106', &
  'I00FW2W202','616BEHW231','I00FW2W242','I00FW2W260','I00FW2W260', &
  'I00FW2W260','200DVEW746','101DVEW375','400DVEW475','101DVEW740', &
  '102DVEW740','101DVEW740','102DVEW740','102DVEW746','102DVEW740', &
  '102DVEW740','200DVEW746','I00FW2W260','203FW2W122','102DVEW740', &
  '101DVEW108','100JB2W202','101DVEW202','102DVEW202','100FW2W202', &
  '103DVEW202','101DVEW017','100JB2W017','101DVEW060','100JB2W108', &
  '104DVEW108','103DVEW108','102DVEW260','100JB2W122','104DVEW122', &
  '102DVEW122','105DVEW122','106DVEW122','101DVEW122','102DVEW019', &
  '101DVEW019','100JB2W073','102DVEW073','101DVEW073','101DVEW999', &
  '101DVEW260','102DVEW101','I12FW2W012','I11FW2W012','I13FW2W012', &
  'I14FW2W012','I21FW2W012','I22FW2W012','102DVEW070','101DVEW240', &
  '102DVEW240','102DVEW119','101DVEW119','103DVEW122','100JB2W019', &
  '101DVEW090'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    DO I=1,81
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
    IF(VOLEQ(1:7).EQ.'616BEHW')THEN
!
!  FOUND VALID WESTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF

    IF(VOLEQ(1:7).EQ.'632BEHW')THEN
!
!  FOUND VALID EASTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF

    IF(VOLEQ(1:7).EQ.'B00BEHW')THEN
!
!  FOUND VALID BLM/BIA/INDUSTRY LANDS EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    RETURN
  ENDIF
!
! Added species 70, 90 and 260 for CruiseProcessing (YW 09/09/2024)
  IF(SPEC.EQ.70) SPEC = 73
  IF(SPEC.EQ.90) SPEC = 93
  IF(SPEC.EQ.260) SPEC = 263

  LAST = 38
  READ(FORST,'(I2)')FORNUM
  DONE = 0

  IF(SPEC.EQ.122 .AND. FORNUM.EQ.8) THEN
     VOLEQ = EQNUM(39)
  ELSEIF(SPEC.EQ.101)THEN
     IF((VAR.EQ.'EM').OR.(VAR.EQ.'em').OR. &
        (VAR.EQ.'IE').OR.(VAR.EQ.'ie').OR. &
       (VAR.EQ.'CI').OR.(VAR.EQ.'ci'))THEN
     VOLEQ = EQNUM(1)
    ENDIF
  ELSEIF((VAR.EQ.'EM'.OR.VAR.EQ.'em').AND.((SPEC.EQ.745).OR. &
         (SPEC.EQ.747).OR.(SPEC.EQ.749)))THEN
     VOLEQ=EQNUM(40)
  ELSE
     FIRST = 1

     DO 5, WHILE (DONE.EQ.0)
        HALF = (LAST - FIRST +1)/2 + FIRST
        IF(FIA(HALF) .EQ. SPEC)THEN
           DONE = HALF
        ELSEIF(FIRST .EQ. LAST) THEN
           ERRFLAG = 1
           DONE = -1
        ELSE IF (FIA(HALF) .LT. SPEC) THEN
           FIRST = HALF
        ELSE
           LAST = HALF - 1
        ENDIF
5 CONTINUE

     IF(DONE .LT. 0)THEN
        IF(SPEC.LT.300)THEN
!              Other Softwood
           VOLEQ = EQNUM(26)
        ELSE
!              Other Hardwood
           VOLEQ = EQNUM(36)
        ENDIF
     ELSE
        VOLEQ = EQNUM(DONE)
     ENDIF
  ENDIF

  RETURN
  END
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R2_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST
  INTEGER SPEC,ERRFLAG
  CHARACTER*10 EQNUM(75)
  INTEGER FIA(42), FIRST, HALF, LAST, DONE ,FORNUM,I
!     SPECIES
!     White fir,        Grand fir,         Corkbark fir,      Subalpine fir,        Juniper,
!     Alligator juniper,Utah Juniper,      Rocky Mtn. Juniper,Eastern redcedar,     Oneseed Juniper,
!     Western Larch,    Engelmann's spruce,White spruce,      Blue spruce,          Whitebark pine,
!     Bristlecode pine, Pinyon Pine,       Lodgepole pine,    Limber pine,          Southwestern white pine,
!     Chihuahua pine,   Ponderosa pine,    Singleleaf pinyon, Border pinyon,        Arizona pinyon,
!     Douglas fir,      Western Redcedar,  Mountain Hemlock,  Other softwoods,      Paper birch,
!     Mountain Mahogany,Cottonwoods,       Plains cottonwood, Quaking aspen,        Narrowleaf cottonwood,
!     Oak,              Arizona white oak, Emory Oak,         Gambel Oak,           Bur Oak,
!     Silverleaf oak,   Other Hardwoods

  DATA (FIA(I),I=1,42)/ 15,  17,  18,  19,  57, &
                        63,  65,  66,  68,  69, &
                        73,  93,  94,  96, 101, &
                       102, 106, 108, 113, 114, &
                       118, 122, 133, 134, 143, &
                       202, 242, 264, 299, 375, &
                       475, 740, 745, 746, 749, &
                       800, 803, 810, 814, 823, &
                       843, 998/
! Changed the default equation for species 66 from 300DVEW060 to 200DVEW066
! and species 814 default from 300DVEW800 to 200DVEW814 to match CruiseProcessing default (YW 9/9/2014)
  DATA (EQNUM(I),I=1,75)/ &
  '200FW2W015','I00FW2W019','I00FW2W019','I00FW2W019','300DVEW060', &
  '300DVEW060','200DVEW065','200DVEW066','300DVEW060','200DVEW069', &
  '407FW2W093','407FW2W093','407FW2W093','407FW2W093','200FW2W122', &
  '200FW2W122','200DVEW106','200FW2W108','200FW2W122','200FW2W122', &
  '200FW2W122','200FW2W122','300DVEW106','300DVEW106','300DVEW106', &
  '200FW2W202','407FW2W093','407FW2W093','300DVEW060','300DVEW999', &
  '200DVEW475','300DVEW999','300DVEW999','200FW2W746','300DVEW999', &
  '300DVEW800','300DVEW800','300DVEW800','200DVEW814','200DVEW823', &
  '300DVEW800','200DVEW998','203FW2W122','213FW2W122','202FW2W108', &
  '200CZ2W202','200DVEW093','200CZ2W093','210DVEW093','200DVEW814', &
  '200DVEW108','210DVEW108','200CZ2W108','212DVEW122','213DVEW122', &
  '200CZ2W122','203CZ2W122','210DVEW122','200DVEW122','203DVEW122', &
  '200CZ2W746','200DVEW746','210DVEW746','200DVEW066','200CZ2W019', &
  '200CZ2W015','200DVEW015','200CZ3W202','200CZ3W093','200CZ3W108', &
  '200CZ3W122','203CZ3W122','200CZ3W746','200CZ3W019','200CZ3W015'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    DO I=1,75
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  DONE = 0

  LAST = 42

  IF(SPEC.EQ.122 .AND. FORNUM.EQ.3) THEN
     VOLEQ = EQNUM(43)
  ELSE IF(SPEC.EQ.122 .AND. FORNUM.EQ.13) THEN
     VOLEQ = EQNUM(44)
  ELSE IF(SPEC.EQ.108 .AND. ((FORNUM.EQ.2).OR.(FORNUM.EQ.14))) THEN
     VOLEQ = EQNUM(45)
  ELSE
     FIRST = 1

     DO 5, WHILE (DONE.EQ.0)
       HALF = (LAST - FIRST +1)/2 + FIRST
        IF(FIA(HALF) .EQ. SPEC)THEN
           DONE = HALF
        ELSEIF(FIRST .EQ. LAST) THEN
           ERRFLAG = 1
           DONE = -1
           ELSE IF (FIA(HALF) .LT. SPEC) THEN
           FIRST = HALF
        ELSE
           LAST = HALF - 1
        ENDIF
5 CONTINUE
     IF(DONE .LT. 0)THEN
!           Other Hardwood
         VOLEQ = EQNUM(42)
     ELSE
        VOLEQ = EQNUM(DONE)
     ENDIF
   ENDIF
!
  RETURN
  END
!
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R3_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST
  INTEGER SPEC,ERRFLAG
  CHARACTER*10 EQNUM(52)
  INTEGER FIA(45), FIRST, HALF, LAST, DONE, FORNUM,I

!     SPECIES
!     White fir,              Grand fir,        Corkbark fir,         Subalpine fir,     Juniper,
!     Juniper,                Alligator juniper,Utah juniper,         Rocky Mtn. Juniper,Eastern redcedar
!     Oneseed juniper,        Western Larch,    Engelmann's spruce,   White spruce,      Blue spruce,
!     Whitebark pine,         Bristlecode pine, Pinyon Pine,          Lodgepole pine,    Limber pine,
!     Southwestern white pine,Chihuahua pine,   Ponderosa pine,       Singleleaf pinyon, Border pinyon,
!     Arizona pinyon,         Douglas fir,      Western Redcedar,     Mountain Hemlock,  Other softwoods,
!     Maple,                  Black maple,      Paper Birch,          Mountain Mahogany, Cottonwoods,
!     Plains cottonwood,      Quaking aspen,    Narrowleaf cottonwood,Oak,               Arizona white oak,
!     Emory oak,              Gambel oak,       Bur oak,              Silverleaf oak,    Other Hardwoods


  DATA (FIA(I),I=1,45)/ 15,  17,  18,  19,  57, &
                        60,  63,  65,  66,  68, &
                        69,  73,  93,  94,  96, &
                       101, 102, 106, 108, 113, &
                       114, 118, 122, 133, 134, &
                       143, 202, 242, 264, 299, &
                       310, 314, 375, 475, 740, &
                       745, 746, 749, 800, 803, &
                       810, 814, 823, 843, 998/
!Updated R3 defaults for the following species (YW 9/9/2024)
  ! Species FROM        TO
  !  15     301DVEW015  NVB0000015
  !  18     300DVEW093  NVB0000015
  !  93     300DVEW093  NVBM330093
  ! 114     300DVEW113  NVBM240119
  ! 202     301DVEW202  300FW2W202
  ! 746     300DVEW746  NVB0000746
  DATA (EQNUM(I),I=1,52)/ &
  'NVB0000015','300DVEW093','NVB0000015','300DVEW093','300DVEW060', &
  '300DVEW060','300DVEW060','300DVEW060','300DVEW060','300DVEW060', &
  '300DVEW060','301DVEW015','NVBM330093','300DVEW093','300DVEW093', &
  '300DVEW113','300DVEW113','300DVEW106','301DVEW202','300DVEW113', &
  'NVBM240119','300DVEW122','300FW2W122','300DVEW106','300DVEW106', &
  '300DVEW106','300FW2W202','301DVEW015','301DVEW015','300DVEW060', &
  '300DVEW310','300DVEW314','300DVEW999','300DVEW999','300DVEW999', &
  '300DVEW999','NVB0000746','300DVEW999','300DVEW800','300DVEW800', &
  '300DVEW800','300DVEW800','300DVEW800','300DVEW800','300DVEW999', &
  '302DVEW202','302DVEW202','302DVEW015','301FW2W122','301FW2W202', &
  '301FW2W015','301FW2W108'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    DO I=1,52
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  DONE = 0

  LAST = 45
! Remove Forest specific default for R3 (YW 9/9/2024)
!      IF(SPEC.EQ.202.AND.(FORNUM.EQ.2.OR.FORNUM.EQ.3.OR.
!     >                    FORNUM.EQ.7)) THEN
!         DONE=46
!      ELSEIF(SPEC.EQ.15.AND.(FORNUM.EQ.2.OR.FORNUM.EQ.3.OR.
!     >                       FORNUM.EQ.7)) THEN
!         DONE=48
!    using Fleweling profile model for Ponderosa pine in following forest:
!    Apache Sitgreaves, Coconino, Kaibab, Tonto. The 300FW2W122 equation
!    is default for PP, so tht other forests 2,3,5,6,8,9,10 and 11 are
!    set to the DVE equation blow.
!      ELSEIF(SPEC.EQ.122.AND.(FORNUM.EQ.2.OR.FORNUM.EQ.3.OR.
!     >                        FORNUM.EQ.5.OR.FORNUM.EQ.6.OR.
!     >                        FORNUM.EQ.8.OR.FORNUM.EQ.9.OR.
!     >                        FORNUM.EQ.11))THEN
!         DONE=22
!    Added profile model for DF, PP, WF and WP to Santa Fe NF
!      ELSEIF(SPEC.EQ.122.AND.FORNUM.EQ.10) THEN
!         DONE = 49
!      ELSEIF(SPEC.EQ.202.AND.FORNUM.EQ.10) THEN
!         DONE = 50
!      ELSEIF(SPEC.EQ.15.AND.FORNUM.EQ.10) THEN
!         DONE = 51
!      ELSEIF(SPEC.EQ.114.AND.FORNUM.EQ.10) THEN
!         DONE = 52
!      ELSE
     FIRST = 1

      DO 5, WHILE (DONE.EQ.0)
        HALF = (LAST - FIRST +1)/2 + FIRST
         IF(FIA(HALF) .EQ. SPEC)THEN
            DONE = HALF
         ELSEIF(FIRST .EQ. LAST) THEN
            ERRFLAG = 1
            DONE = -1
            ELSE IF (FIA(HALF) .LT. SPEC) THEN
            FIRST = HALF
         ELSE
            LAST = HALF - 1
         ENDIF
5 CONTINUE
!      ENDIF
  IF(DONE .LT. 0)THEN
!        Unknown
     VOLEQ = EQNUM(45)
  ELSE
     VOLEQ = EQNUM(DONE)
  ENDIF
!
  RETURN
  END
!
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R4_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST
  INTEGER SPEC,ERRFLAG
  CHARACTER*10 EQNUM(57)
  INTEGER FIA(27), FIRST, HALF, LAST, DONE, FORNUM,I

! SPECIES
! California red fir,    Juniper,          Western juniper,      Rocky Mtn juniper,    Western larch,
! Incense-cedar,         Common pinyon,    Sugar pine,           Western white pine,   Singleleaf pinyon,
! Pacific yew,           Western redcedar, Western hemlock,      Mountain hemlock,     Other softwoods,
! Box elder,             Rocky Mtn maple,  Bigtooth maple,       Curlleaf Mtn Mahogany,Quaking aspen,
! Black cottonwood,      Fremont cottonwod,Narrowleaf cottonwood,Oak (sp.)             Gambel oak,
! Other hardwoods,       Other

  DATA (FIA(I),I=1,27)/   20,   60,   64,   66,   73, &
                          81,  106,  117,  119,  133, &
                         231,  242,  263,  264,  299, &
                         313,  321,  322,  475,  746, &
                         747,  748,  749,  800,  814, &
                         998,  999/

  DATA (EQNUM(I),I=1,57)/ &
  '400MATW020','300DVEW060','300DVEW060','400DVEW066','400MATW073', &
  '400MATW081','300DVEW106','400MATW117','400MATW117','400DVEW133', &
  '400DVEW998','400MATW081','400MATW015','401MATW015','400MATW108', &
  '400MATW108','400MATW108','400MATW108','400DVEW475','400MATW746', &
  '400DVEW998','400MATW108','400MATW108','300DVEW800','300DVEW800', &
  '400DVEW998','400DVEW998','I15FW2W017','401MATW015','400MATW015', &
  'I15FW2W017','405MATW019','400MATW019','401DVEW065','400DVEW065', &
  'I15FW2W093','407FW2W093','400MATW093','401MATW108','400MATW108', &
  'I15FW2W122','401MATW122','402MATW122','403MATW122','400MATW122', &
  'I15FW2W202','405MATW202','400MATW202','400DVEW064','400DVEW106', &
  '407MATW093','401MATW202','407MATW093','401MATW202','401MATW202', &
  '407MATW093','407MATW093'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    DO I=1,57
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  DONE = 0

  LAST = 27

!     White fir
  IF(SPEC.EQ.15) THEN
     IF (FORNUM.EQ.2.OR.FORNUM.EQ.12.OR.FORNUM.EQ.13.OR. &
                                                 FORNUM.EQ.6)THEN
        VOLEQ = EQNUM(28)
     ELSEIF (FORNUM.EQ.9.OR.FORNUM.EQ.17)THEN
        VOLEQ = EQNUM(29)
     ELSE
        VOLEQ = EQNUM(30)
     ENDIF
!     Grand fir
  ELSEIF(SPEC.EQ.17) THEN
     IF (FORNUM.EQ.2.OR.FORNUM.EQ.12.OR.FORNUM.EQ.13)THEN
        VOLEQ = EQNUM(31)
     ELSE
        VOLEQ = EQNUM(30)
     ENDIF
!     Subalpine fir
  ELSEIF(SPEC.EQ.19)THEN
     IF (FORNUM.EQ.5)THEN
        VOLEQ = EQNUM(32)
     ELSE
        VOLEQ = EQNUM(33)
     ENDIF
!     Utah Juniper
  ELSEIF(SPEC.EQ.64) THEN
     IF((FORNUM.EQ.1).OR.(FORNUM.EQ.4).OR.(FORNUM.EQ.7).OR. &
        (FORNUM.EQ.8).OR.(FORNUM.EQ.9).OR.(FORNUM.EQ.10).OR. &
        (FORNUM.EQ.17).OR.(FORNUM.EQ.18).OR.(FORNUM.EQ.19))THEN
        VOLEQ=EQNUM(49)
     ELSEIF((FORNUM.EQ.2).OR.(FORNUM.EQ.6).OR.(FORNUM.EQ.12).OR. &
        (FORNUM.EQ.13).OR.(FORNUM.EQ.14))THEN
        VOLEQ=EQNUM(49)
     ENDIF
  ELSEIF(SPEC.EQ.65) THEN
     IF (FORNUM.EQ.3.OR.FORNUM.EQ.5.OR. &
         FORNUM.EQ.15.OR.FORNUM.EQ.16) THEN
        VOLEQ = EQNUM(34)
      ELSE
        VOLEQ = EQNUM(35)
      ENDIF
!     Engelmann spruce
!     Blue spruce
  ELSEIF(SPEC.EQ.93.OR.SPEC.EQ.96) THEN
     IF (FORNUM.EQ.2.OR.FORNUM.EQ.12.OR.FORNUM.EQ.13)THEN
        VOLEQ = EQNUM(36)
     ELSEIF (FORNUM.EQ.7) THEN
        VOLEQ = EQNUM(37)
     ELSEIF (FORNUM.EQ.8) THEN
        VOLEQ = EQNUM(56)
      ELSE
        VOLEQ = EQNUM(38)
      ENDIF
!     Whitebark pine
!     Lodgepole pine
!     Limber pine
!     Bristlecone pine
  ELSEIF(SPEC.EQ.101.OR.SPEC.EQ.108.OR.SPEC.EQ.113.OR. &
         SPEC.EQ.142) THEN
     IF (FORNUM.EQ.9.OR.FORNUM.EQ.17) THEN
        VOLEQ = EQNUM(39)
      ELSE
        VOLEQ = EQNUM(40)
      ENDIF
!     pinyon pine
  ELSEIF(SPEC.EQ.106) THEN
     IF((FORNUM.EQ.1).OR.(FORNUM.EQ.4).OR.(FORNUM.EQ.7).OR. &
        (FORNUM.EQ.8).OR.(FORNUM.EQ.9).OR.(FORNUM.EQ.10).OR. &
        (FORNUM.EQ.17).OR.(FORNUM.EQ.18).OR.(FORNUM.EQ.19))THEN
        VOLEQ=EQNUM(50)
     ENDIF
!     Ponderosa pine
  ELSEIF(SPEC.EQ.122) THEN
     IF (FORNUM.EQ.2.OR.FORNUM.EQ.12.OR.FORNUM.EQ.13)THEN
        VOLEQ = EQNUM(41)
     ELSEIF (FORNUM.EQ.01)THEN
        VOLEQ = EQNUM(42)
     ELSEIF (FORNUM.EQ.7.OR.FORNUM.EQ.8.OR.FORNUM.EQ.10.OR. &
             FORNUM.EQ.18.OR.FORNUM.EQ.19) THEN
        VOLEQ = EQNUM(43)
     ELSEIF (FORNUM.EQ.9.OR.FORNUM.EQ.17)THEN
        VOLEQ = EQNUM(44)
      ELSE
        VOLEQ = EQNUM(45)
      ENDIF
!     Douglas fir
  ELSEIF(SPEC.EQ.202) THEN
     IF (FORNUM.EQ.2.OR.FORNUM.EQ.12.OR.FORNUM.EQ.13)THEN
        VOLEQ = EQNUM(46)
     ELSEIF (FORNUM.EQ.05)THEN
        VOLEQ = EQNUM(47)
      ELSE
        VOLEQ = EQNUM(48)
      ENDIF
!     Other hardwoods
  ELSEIF(SPEC.EQ.998) THEN
    IF((FORNUM.EQ.2).OR.(FORNUM.EQ.6).OR.(FORNUM.EQ.12).OR. &
          (FORNUM.EQ.13).OR.(FORNUM.EQ.14))THEN
          VOLEQ=EQNUM(17)
    ELSE
          VOLEQ=EQNUM(21)
    ENDIF
  ELSE
     FIRST = 1

      DO 5, WHILE (DONE.EQ.0)
        HALF = (LAST - FIRST +1)/2 + FIRST
         IF(FIA(HALF) .EQ. SPEC)THEN
            DONE = HALF
         ELSEIF(FIRST .EQ. LAST) THEN
            ERRFLAG = 1
            DONE = -1
            ELSE IF (FIA(HALF) .LT. SPEC) THEN
            FIRST = HALF
         ELSE
            LAST = HALF - 1
         ENDIF
5 CONTINUE

     IF(DONE .LT. 0)THEN
!           Other Hardwood
        VOLEQ = EQNUM(26)
     ELSE
        VOLEQ = EQNUM(DONE)
     ENDIF
   ENDIF

  RETURN
  END


!//////////////////////////////////////////////////////////////////
  SUBROUTINE R5_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST,VAR
  INTEGER SPEC,ERRFLAG
  CHARACTER*10 EQNUM(108)
  INTEGER FIA(72), FIRST, HALF, LAST, DONE, FORNUM,I

!     SPECIES
!     Pacific silver fir,    White fir,             Grand fir,           Subalpine fir,       California red fir,
!     Shasta red fir,        Noble fir,             Port Orford cedar,   California juniper,  Juniper,
!     Utah juniper,          Western Larch,         Incense cedar,       Brewer spruce,       Engelmann spruce,
!     White bark pine,       Knobcone pine,         Foxtail pine,        Lodgepole pine,      Coulter pine,
!     Limber pine,           Jeffrey pine,          Sugar pine,          Western white pine,  Ponderosa pine,
!     Monterey pine,         Grey pine,             Singleleaf pinyon,   Washoe pine,         Great Basin bristlecone pine
!     Bigcone Douglas-fir,   Douglas-fir,
!     Redwood,               Giant sequoia,         Pacific yew,         Western red cedar,   California nutmeg,
!     Western hemlock,       Mountain hemlock,      Other softwoods,     Koa,                 Bigleaf maple,
!     California buckeye,    Red alder,             White alder          Pacific madrone,     Golden chinkapin,
!     Curl-leaf Mtn Mahogany,Birchleaf Mtn Mahogany,Pacific dogwood,     Eucalyptus,          Robusta Eucaltyptus,    Oregon Ash,
!     Walnut,                Tanoak,                Ohia,                California sycamore, Quacking aspen,
!     Black cottonwood,      Bitter Cherry,         California live oak, Canyon live oak,     Blue oak,
!     Engelmann's oak,       Oregon white oak,      California black oak,California white oak,Interior live oak,
!     Willow,                California laurel,     Other hardwoods,     unkown
! Changed the default for species 17 to be same as CruiseProcessing default and added species 514
! Species From       To
!  17     500WO2W015 I15FW2W017
! 514                H00SN2W510
  DATA (FIA(I),I=1,32)/   11,   15,   17,   19,   20, &
                          21,   22,   41,   62,   64, &
                          65,   73,   81,   92,   93, &
                         101,  103,  104,  108,  109, &
                         113,  116,  117,  119,  122, &
                         124,  127,  133,  137,  142, &
                         201,  202/
  DATA (FIA(I),I=33,72)/ 211,  212,  231,  242,  251, &
                         263,  264,  299,  301,  312, &
                         333,  351,  352,  361,  431, &
                         475,  478,  492,  510,  514, 542, &
                         600,  631,  671,  730,  746, &
                         747,  768,  801,  805,  807, &
                         811,  815,  818,  821,  839, &
                         920,  981,  998,  999/

  DATA (EQNUM(I),I=1,108)/ &
  '500WO2W015','500WO2W015','I15FW2W017','500WO2W020','500WO2W020', &
  '500WO2W020','500WO2W020','500WO2W081','500DVEW060','500DVEW060', &
  '500DVEW060','500WO2W202','500WO2W081','500WO2W015','500WO2W015', &
  '500WO2W108','500WO2W108','500WO2W108','500WO2W108','500WO2W108', &
  '500WO2W108','500WO2W116','500WO2W117','500WO2W117','500WO2W122', &
  '500WO2W108','500WO2W108','500WO2W116','500WO2W117','500WO2W108', &
  '500WO2W202','500WO2W202', &
  '500WO2W211','500DVEW212','500WO2W108','500WO2W081','500DVEW807', &
  '500WO2W015','500WO2W015','500WO2W108','H00SN2W301','500DVEW312', &
  '500DVEW807','500DVEW351','500DVEW351','500DVEW361','500DVEW431', &
  '500DVEW801','500DVEW801','500DVEW807','H01SN2W510','H00SN2W510', &
   '500DVEW807', &
  '500DVEW818','500DVEW631','H00SN2W671','500DVEW818','500DVEW818', &
  '500DVEW818','500DVEW801','500DVEW801','500DVEW805','500DVEW807', &
  '500DVEW811','500DVEW815','500DVEW818','500DVEW821','500DVEW839', &
  '500DVEW807','500DVEW981','500DVEW981','500DVEW631','532WO2W015', &
  '532WO2W020','532WO2W081','532WO2W108','532WO2W117','532WO2W122', &
  '532WO2W202','H00SN2W510','532WO2W020','532WO2W202','H00SN2W510', &
  '532WO2W081','532WO2W116','616BEHW000','532WO2W108','532WO2W122', &
  '532WO2W211','532WO2W117','532WO2W015','516TRFW015','532TRFW015', &
  '516TRFW021','532TRFW021','516TRFW060','532TRFW060','516TRFW081', &
  '532TRFW081','516TRFW108','532TRFW108','516TRFW117','532TRFW117', &
  '516TRFW122','532TRFW122','516TRFW202','616TRFW211','632TRFW211'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    DO I=1,107
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  DONE = 0
  IF(SPEC.EQ.60) SPEC = 62
!     Whitebark pine
  IF(SPEC.EQ.101) THEN
      IF(VAR.EQ.'SO' .OR. VAR.EQ.'so') THEN
           DONE=30
      ELSE
           DONE=40
      ENDIF
!     other softwoods
  ELSEIF(SPEC.EQ.299 .OR. SPEC.EQ.290) THEN
      IF(VAR.EQ.'SO' .OR. VAR.EQ.'so') THEN
          DONE=30
      ELSEIF(VAR.EQ.'NC' .OR. VAR.EQ.'nc') THEN
          DONE=30
      ELSE
          DONE=40
      ENDIF
!     other hardwoods
  ELSEIF(SPEC.EQ.998.) THEN
      IF(VAR.EQ.'SO' .OR. VAR.EQ.'so') THEN
          DONE=70
      ELSEIF(VAR.EQ.'WS' .OR. VAR.EQ.'ws') THEN
          DONE=67
      ELSEIF(VAR.EQ.'NC' .OR. VAR.EQ.'nc') THEN
          DONE=70
      ELSE
          DONE=60
      ENDIF
!      ELSE
   ENDIF
   LAST = 72

   FIRST = 1

   DO 5, WHILE (DONE.EQ.0)
       HALF = (LAST - FIRST +1)/2 + FIRST
       IF(FIA(HALF) .EQ. SPEC)THEN
           DONE = HALF
       ELSEIF(FIRST .EQ. LAST) THEN
           ERRFLAG = 1
           DONE = -1
       ELSE IF (FIA(HALF) .LT. SPEC) THEN
           FIRST = HALF
       ELSE
           LAST = HALF - 1
       ENDIF
5 CONTINUE

   IF(DONE .LT. 0) THEN
       VOLEQ = EQNUM(72)
   ELSE
       VOLEQ = EQNUM(DONE)
   ENDIF
!      ENDIF

  RETURN
  END


!//////////////////////////////////////////////////////////////////
  SUBROUTINE R6_EQN(VAR,FORST,DIST,SPEC,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST,VAR,DIST
  CHARACTER*3 ASPEC
  INTEGER SPEC,ERRFLAG,FORNUM,I,FIA(53),DISTNUM
  CHARACTER*10 EQNUMI(87),EQNUM(44),EQNUMRD(2)
  CHARACTER*10 EQNUMC(39),EQNUMF(50),EQNUMD(15)
  INTEGER FIRST, HALF, LAST, DONEI, DONEF, DONERD

!     All species 32' logs, All species 16' logs,
!     Balsam fir,      Grand Fir,  Alpine Fir,       Western Larch, Engelmann spruce
!     Lodgepole pine,  White Pine, Ponderosa pine,   Douglas Fir,   Western red cedar,
!     Mountain Hemlock,Douglas fir,Western red cedar,Mountain Hemlock

  DATA (FIA(I),I=1,53)/  011, 015, 017, 019, 020, &
                         021, 022, 042, 064, 066, &
                         072, 073, 081, 093, 098, &
                         101, 103, 106, 108, 113, &
                         116, 117, 119, 122, 202, &
                         211, 231, 242, 263, 264, &
                         290, 299, 312, 321, 351, &
                         352, 361, 375, 431, 475, &
                         478, 492, 500, 631, 740, &
                         746, 747, 768, 815, 818, &
                         920, 998, 999/
!
  DATA (EQNUMI(I),I=1,87)/ &
  'I00FW2W012','I00FW2W017','I00FW2W019','I00FW2W073','I00FW2W093', &
  'I00FW2W108','I00FW2W119','I00FW2W122','I00FW2W202','I00FW2W242', &
  'I00FW2W260','I00FW2W263','I11FW2W012','I11FW2W017','I11FW2W019', &
  'I11FW2W073','I11FW2W093','I11FW2W108','I11FW2W119','I11FW2W122', &
  'I11FW2W202','I11FW2W242','I11FW2W260','I11FW2W263','I12FW2W012', &
  'I12FW2W017','I12FW2W019','I12FW2W073','I12FW2W093','I12FW2W108', &
  'I12FW2W119','I12FW2W122','I12FW2W202','I12FW2W242','I12FW2W260', &
  'I12FW2W263','I13FW2W012','I13FW2W017','I13FW2W019','I13FW2W073', &
  'I13FW2W093','I13FW2W108','I13FW2W119','I13FW2W122','I13FW2W202', &
  'I13FW2W242','I13FW2W260','I13FW2W263','I14FW2W012','I14FW2W017', &
  'I14FW2W019','I14FW2W073','I14FW2W093','I14FW2W108','I14FW2W119', &
  'I14FW2W122','I14FW2W202','I14FW2W242','I14FW2W260','I14FW2W263', &
  'I00FW2W012','I00FW2W202','I00FW2W093','I00FW2W017','I00FW2W108', &
  'I00FW2W122','I00FW2W019','I00FW2W073','I00FW2W119','I00FW2W242', &
  'I00FW2W260','I13FW2W202','I11FW2W202','I12FW2W202','I14FW2W017', &
  'I13FW2W017','I14FW2W122','I12FW2W122','I13FW2W122','I11FW2W122', &
  'I11FW2W017','I11FW2W073','I11FW2W242','I12FW2W017','I12FW2W073', &
  'I13FW2W073','I14FW2W073'/

  DATA (EQNUMF(I),I=1,50)/ &
  'F00FW2W202','F00FW2W242','F00FW2W263','F01FW2W202','F01FW2W242', &
  'F01FW2W263','F02FW2W202','F02FW2W242','F02FW2W263','F03FW2W202', &
  'F03FW2W242','F03FW2W263','F04FW2W202','F04FW2W242','F04FW2W263', &
  'F05FW2W202','F05FW2W242','F05FW2W263','F06FW2W202','F06FW2W242', &
  'F06FW2W263','F07FW2W202','F07FW2W242','F07FW2W263','F08FW2W202', &
  'F08FW2W242','F08FW2W263','F03FW2W202','F01FW2W202','F02FW2W202', &
  'F00FW2W202','F04FW2W202','F08FW2W202','F07FW2W202','F06FW2W202', &
  'F05FW2W202','F00FW2W242','F00FW2W260','F01FW2W242','F01FW2W260', &
  'F02FW2W242','F02FW2W260','F03FW2W242','F03FW2W260','F04FW2W242', &
  'F04FW2W260','F05FW2W260','F06FW2W260','F07FW2W260','F08FW2W260'/

  DATA (EQNUMC(I),I=1,39)/ &
  'I21FW2W012','I21FW2W017','I21FW2W019','I21FW2W073','I21FW2W093', &
  'I21FW2W108','I21FW2W119','I21FW2W122','I21FW2W202','I21FW2W242', &
  'I21FW2W260','I22FW2W012','I22FW2W017','I22FW2W019','I22FW2W073', &
  'I22FW2W093','I22FW2W108','I22FW2W119','I22FW2W122','I22FW2W202', &
  'I22FW2W242','I22FW2W260','I23FW2W012','I23FW2W017','I23FW2W019', &
  'I23FW2W073','I23FW2W093','I23FW2W108','I23FW2W119','I23FW2W122', &
  'I23FW2W202','I23FW2W242','I23FW2W260','I22FW2W122','I21FW2W122', &
  'I23FW2W122','I21FW2W073','I22FW2W073','I23FW2W073'/

  DATA (EQNUMD(I),I=1,15)/ &
  '601DVEW205','601DVEW263','601DVEW015','602DVEW122','602DVEW204', &
  '602DVEW015','602DVEW108','602DVEW122','601DVEW205','602DVEW204', &
  '602DVEW108','600DVEW122','601DVEW015','601DVEW263','602DVEW015'/

  DATA (EQNUM(I),I=1,44)/ &
  'B00BEHW011','B00BEHW015','B00BEHW017','B00BEHW015','B00BEHW021', &
  'B00BEHW021','B00BEHW022','B00BEHW041','B00BEHW042','B00BEHW242', &
  'B00BEHW073','B00BEHW073','B00BEHW081','B00BEHW093','B00BEHW098', &
  'B00BEHW119','B00BEHW108','B00BEHW108','B00BEHW116','B00BEHW117', &
  'B00BEHW119','B00BEHW122','B01BEHW202','B00BEHW211','B00BEHW231', &
  'B00BEHW242','B00BEHW260','B00BEHW260','B00BEHW260','B00BEHW312', &
  'B00BEHW351','B00BEHW361','B00BEHW361','B00BEHW431','B00BEHW542', &
  'B00BEHW631','B00BEHW747','B00BEHW800','B00BEHW800','B00BEHW998', &
  'B00BEHW999','B02BEHW202','B03BEHW202','B01BEHW202'/

  !Equation for Red Alder in GP and Siuslaw
  DATA (EQNUMRD(I),I=1,2)/ &
  'A16CURW351','NVBM240351'/


!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    IF(VOLEQ(1:7).EQ.'616BEHW')THEN
!
!  FOUND VALID WESTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF

    if(VOLEQ(1:7).EQ.'628BEHW')THEN
        SPEC=8888
        RETURN
    ENDIF

    IF(VOLEQ(1:7).EQ.'632BEHW')THEN
!
!  FOUND VALID EASTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
!
!  CHECK FOR VALID PNW TARRIF EQUATION
!
    IF(VOLEQ(1:7).EQ.'616TRFW')THEN
      SPEC=8888
      RETURN
    ENDIF

    IF(VOLEQ(1:7).EQ.'632TRFW')THEN
     SPEC=8888
     RETURN
   ENDIF

    DO I=1,44
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID BLM/BIA/INDUSTRY LANDS EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,87
    IF(VOLEQ.EQ.EQNUMI(I))THEN
!
!  FOUND VALID INGY EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,50
    IF(VOLEQ.EQ.EQNUMF(I))THEN
!
!  FOUND VALID WESTSIDE FLEWELLING EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,39
    IF(VOLEQ.EQ.EQNUMC(I))THEN
!
!  FOUND VALID INGY EQUATION NUMBER - CANADIAN MODEL
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,15
    IF(VOLEQ.EQ.EQNUMD(I))THEN
!
!  FOUND VALID DIRECT VOLUME ESTIMATORS
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
!2025/11/24 added to check Red Alder EQNUMRD
    DO I=1,2
        IF(VOLEQ.EQ.EQNUMRD(I))THEN
            ! FOUND EQUATION FOR RED ALDER
            SPEC=8888
            RETURN
        ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  READ(DIST,'(I2)')DISTNUM
  DONEI = 0
  DONEF = 0
  DONERD = 0

!     Westside Variants
  IF(VAR.EQ.'PN' .OR. VAR.EQ.'WC' .OR. VAR.EQ.'NC' .OR. &
     VAR.EQ.'CA' .OR. VAR.EQ.'OC' .OR. VAR.EQ.'OP')THEN

!        Gifford Pinchot
     IF(FORNUM.EQ.3)THEN
        IF(SPEC.EQ. 11) THEN
            DONEI = 26
        ELSE IF (SPEC.EQ.19) THEN
            DONEI=6
        ELSE IF(SPEC.EQ.263 .OR. SPEC.EQ.260) THEN
            DONEF = 3
            IF(DISTNUM.EQ.1) THEN
              DONEF = 21
            ELSE IF (DISTNUM.EQ.5) THEN
               DONEF = 15
            ENDIF
        ELSE IF(SPEC.EQ.202) THEN
            DONEF = 10
            IF(DISTNUM.EQ.1) THEN
              DONEF = 22
            ELSE IF (DISTNUM.EQ.5) THEN
               DONEF = 10
            ENDIF
        ELSE IF(SPEC.EQ.351) THEN
            DONERD = 1
        ENDIF
!        Mt Hood
     ELSE IF(FORNUM.EQ.6)THEN
        IF(SPEC.EQ.11) THEN
            DONEI = 26
         ELSE IF(SPEC.EQ.17) THEN
            DONEI = 38
         ELSE IF(SPEC.EQ.93)THEN
            DONEI = 17
         ELSE IF(SPEC.EQ.108)THEN
            DONEI = 18
         ELSE IF(SPEC.EQ.122)THEN
            DONEI = 32
         ELSE IF(SPEC.EQ.263 .OR. SPEC.EQ.260) THEN
            DONEI = 23
         ELSE IF(SPEC.EQ.22) THEN
            DONEI = 38
         ELSE IF(SPEC.EQ.202) THEN
            DONEF = 10
            IF(DISTNUM.EQ.5 .OR. DISTNUM.EQ.9) DONEF = 10
          END IF
!        Mt Baker - Snoqualmie
     ELSE IF(FORNUM.EQ.5)THEN
         IF(SPEC.EQ.263 .OR. SPEC.EQ.260) THEN
            DONEF = 12
         ELSE IF(SPEC.EQ.202) THEN
            DONEF = 25
         ENDIF
!        Rogue River/Siskiyo
     ELSE IF(FORNUM.EQ.10.OR.FORNUM.EQ.11) THEN
         IF(SPEC.EQ.15) THEN
            IF((VAR.NE.'NC').OR.(VAR.NE.'nc'))THEN
               DONEI = 5
            ENDIF
         ELSE IF(SPEC.EQ.122) THEN
            IF((VAR.NE.'NC').OR.(VAR.NE.'nc'))THEN
               DONEI = 4
            ENDIF
         ELSE IF(SPEC.EQ.202) THEN
            DONEF = 19
         ELSE IF(SPEC.EQ.263) THEN
            DONEF = 21
         ENDIF
!        Siuslaw
     ELSE IF(FORNUM.EQ.12) THEN
         IF(SPEC.EQ.202) THEN
            DONEF = 1
         ELSE IF(SPEC.EQ.263) THEN
            DONEF = 12
         ELSE IF(SPEC.EQ.351) THEN
            DONERD = 2
        ENDIF
!        Olympic
     ELSE IF(FORNUM.EQ.9) THEN
         IF(SPEC.EQ.202) THEN
            DONEF = 10
         ELSE IF(SPEC.EQ.98) THEN
            DONEF = 12
         ELSE IF(SPEC.EQ.263) THEN
            DONEF = 3
         ENDIF
!        Umpqua
     ELSE IF(FORNUM.EQ.15) THEN
         IF(SPEC.EQ.15) THEN
            DONEI = 2
         ELSE IF(SPEC.EQ.20)THEN
            DONEI = 1
         ELSE IF(SPEC.EQ.81)THEN
            DONEI = 4
         ELSE IF(SPEC.EQ.93)THEN
            DONEI = 5
         ELSE IF(SPEC.EQ.108)THEN
            DONEI = 6
         ELSE IF(SPEC.EQ.122)THEN
            DONEI = 4
         ELSE IF(SPEC.EQ.202)THEN
            DONEF = 1
         ELSE IF(SPEC.EQ.242.OR.SPEC.EQ.21)THEN
            DONEI = 1
         ELSE IF(SPEC.EQ.263)THEN
            DONEI = 23
         ELSE IF(SPEC.EQ.264)THEN
            DONEI = 10
         ELSE IF(SPEC.EQ.11)THEN
            DONEI = 76
         ELSE IF(SPEC.EQ.103)THEN
            DONEI = 11
         ENDIF
!         Willamette
     ELSE IF(FORNUM.EQ.18) THEN
         IF(SPEC.EQ.22) THEN
            DONEI = 6
         ELSE IF(SPEC.EQ.202)THEN
            DONEF = 16
         ELSE IF(SPEC.EQ.263)THEN
            DONEF = 12
        !Added Incense Cedar (81) to use I00FW2W073 and Grand fir (17) to use I00FW2W017 (20250227)
         ELSE IF(SPEC.EQ.17)THEN
            DONEI = 64
         ELSE IF(SPEC.EQ.81)THEN
            DONEI = 68
         ENDIF
     ENDIF

     IF(DONEI.GT.0) THEN
         VOLEQ = EQNUMI(DONEI)
      ELSE IF(DONEF.GT.0) THEN
         VOLEQ = EQNUMF(DONEF)
      ELSE IF(DONERD.GT.0) THEN
          VOLEQ = EQNUMRD(DONERD)
      ELSE
!           No INGY, find Behre's hyperbola model
        LAST = 53

        FIRST = 1

        DO 5, WHILE (DONEF.EQ.0)
           HALF = (LAST - FIRST +1)/2 + FIRST
           IF(FIA(HALF) .EQ. SPEC)THEN
              DONEF = HALF
           ELSEIF(FIRST .EQ. LAST) THEN
              ERRFLAG = 1
              DONEF = -1
           ELSE IF (FIA(HALF) .LT. SPEC) THEN
              FIRST = HALF
           ELSE
              LAST = HALF - 1
           ENDIF
5 CONTINUE

        IF(DONEF .LT. 0) THEN
           VOLEQ = '616BEHW000'
        ELSE
            IF(FIA(DONEF) .LT. 10) THEN
              WRITE(ASPEC,'(I1)')FIA(DONEF)
              VOLEQ(1:9) = '616BEHW00'
              VOLEQ(10:10) = ASPEC
            ELSE IF (FIA(DONEF) .LT. 100) THEN
              WRITE(ASPEC,'(I2)')FIA(DONEF)
              VOLEQ(1:8) = '616BEHW0'
              VOLEQ(9:10) = ASPEC
            ELSE
              WRITE(ASPEC,'(I3)')FIA(DONEF)
              VOLEQ(1:7) = '616BEHW'
              VOLEQ(8:10) = ASPEC
            ENDIF
        ENDIF
      ENDIF

!     Eastside Variants
  ELSE
!        first check for INGY equations by forest and species
!        Deschutes
     IF(FORNUM.EQ.1) THEN
        IF(SPEC.EQ.11 .OR. SPEC.EQ.15 .OR. SPEC.EQ.17 &
                      .OR. SPEC.EQ.21) THEN
           IF((VAR.NE.'SO').OR.(VAR.NE.'so'))DONEI = 2
        ELSE IF(SPEC.EQ. 73) THEN
           DONEI = 16
        ELSE IF(SPEC.EQ. 108) THEN
           DONEI = 18
        ELSE IF(SPEC.EQ. 122) THEN
           DONEI = 20
        ELSE IF(SPEC.EQ.202 ) THEN
           DONEI = 21
        ELSE IF(SPEC.EQ. 81) THEN
           DONEI = 22
        ELSE IF(SPEC.EQ. 117) THEN
           DONEI = 44
        ENDIF
!        Fremont
     ELSEIF(FORNUM.EQ.2 .OR. FORNUM.EQ.20) THEN
         IF(SPEC.EQ. 15 .OR. SPEC.EQ.17) THEN
           DONEI = 14
         ELSE IF(SPEC.EQ.81) THEN
           DONEI = 9
         ELSE IF(SPEC.EQ.108) THEN
           DONEI = 6
         ELSE IF(SPEC.EQ.122) THEN
           DONEI = 8
         ELSE IF(SPEC.EQ.202) THEN
          DONEI = 2
         ENDIF
!        Gifford Pinchot
     ELSE IF(FORNUM.EQ.3)THEN
        IF(SPEC.EQ. 11) THEN
            DONEI = 26
         ELSE IF(SPEC.EQ.263 .OR. SPEC.EQ.260) THEN
            IF(DISTNUM.EQ.3) DONEF = 3
        ELSE IF(SPEC.EQ.202) THEN
            IF(DISTNUM.EQ.3) DONEI = 10
        ENDIF
!        Malheur
     ELSE IF(FORNUM.EQ.4)THEN
         IF(SPEC.EQ. 108) THEN
           DONEI = 30
         ELSE IF(SPEC.EQ.17 .OR. SPEC.EQ.15) THEN
           DONEI = 26
         ELSE IF(SPEC.EQ.122) THEN
           DONEI = 32
         ELSE IF(SPEC.EQ.202) THEN
           DONEI = 33
         ENDIF
!        Mt Hood
     ELSE IF(FORNUM.EQ.6)THEN
        IF(SPEC.EQ.11) THEN
            DONEI = 26
        ELSE IF(SPEC.EQ.17) THEN
            DONEI = 38
        ELSE IF(SPEC.EQ.93) THEN
            DONEI = 17
        ELSE IF(SPEC.EQ.108) THEN
            DONEI = 18
        ELSE IF(SPEC.EQ.122) THEN
            DONEI = 32
        ELSE IF(SPEC.EQ.263 .OR. SPEC.EQ.260) THEN
            DONEI = 23
        ELSE IF(SPEC.EQ.22) THEN
            DONEI = 38
        ELSE IF(SPEC.EQ.202) THEN
            DONEF=16
            IF(DISTNUM.EQ.1 .OR. DISTNUM.EQ.6)DONEF = 22
          END IF
!        Ochoco
     ELSE IF(FORNUM.EQ.7) THEN
        IF(SPEC.EQ. 17 .OR. SPEC.EQ. 15) THEN
           DONEI = 26
        ELSE IF(SPEC.EQ. 73) THEN
           DONEI = 28
        ELSE IF(SPEC.EQ.122) THEN
           DONEI = 32
        ELSE IF(SPEC.EQ.202) THEN
           DONEI = 33
         ELSE IF(SPEC.EQ.108) THEN
            DONEI = 30
        ENDIF
!        Umatilla
     ELSE IF(FORNUM.EQ.14)THEN
        IF(SPEC.EQ. 17 .OR. SPEC.EQ. 15) THEN
           DONEI = 38
        ELSE IF(SPEC.EQ. 19)THEN
           DONEI = 3
        ELSE IF(SPEC.EQ. 73) THEN
           DONEI = 45
        ELSE IF(SPEC.EQ. 93) THEN
           DONEI = 5
        ELSE IF(SPEC.EQ.108) THEN
           DONEI = 6
        ELSE IF(SPEC.EQ.122) THEN
           DONEI = 44
        ELSE IF(SPEC.EQ.202) THEN
           DONEI = 38
        ENDIF
!        Wallowa-Whitman
     ELSE IF(FORNUM.EQ.16) THEN
        IF(SPEC.EQ. 17 .OR. SPEC.EQ. 15) THEN
           DONEI = 14
        ELSE IF(SPEC.EQ.73) THEN
            DONEI = 40
        ELSE IF(SPEC.EQ.93) THEN
            DONEI = 17
        ELSE IF(SPEC.EQ.108) THEN
            DONEI = 6
        ELSE IF(SPEC.EQ.122) THEN
           DONEI = 20
        ELSE IF(SPEC.EQ.202) THEN
           DONEI = 21
        ENDIF
!        Okanogan - Wenatchee
     ELSE IF(FORNUM.EQ.8 .OR. FORNUM.EQ.17) THEN
        IF(SPEC.EQ. 17) THEN
            DONEI = 14
            IF(DISTNUM.EQ.2 .OR. DISTNUM.EQ.3 .OR. DISTNUM.EQ.5 &
                                    .OR. DISTNUM.EQ.7 ) DONEI = 14
        ELSE IF(SPEC.EQ.202) THEN
            DONEI = 33
            IF(DISTNUM.EQ.2 .OR. DISTNUM.EQ.3 .OR. DISTNUM.EQ.4 &
                .OR. DISTNUM.EQ.5 .OR. DISTNUM.EQ.7 .OR. &
                DISTNUM.EQ.9) DONEI = 33
        ELSE IF(SPEC.EQ.108) THEN
            DONEI = 30
            IF(DISTNUM.EQ.4) DONEI = 30
        ELSE IF(SPEC.EQ.93) THEN
            DONEI = 17
            IF(DISTNUM.EQ.9) DONEI = 17
        ELSE IF(SPEC.EQ.122.OR.SPEC.EQ.73) THEN
            DONEI = 32
            IF(DISTNUM.EQ.4) THEN
              DONEI = 32
            ELSE IF(DISTNUM.EQ.2 .OR. DISTNUM.EQ.3 .OR. &
              DISTNUM.EQ.5.OR. DISTNUM.EQ.7) THEN
              DONEI = 20
            ENDIF
        ENDIF
!        Colville
     ELSE IF(FORNUM.EQ.21) THEN
         IF(SPEC.EQ.17.OR.SPEC.EQ.260) THEN
            DONEI = 14
         ELSE IF(SPEC.EQ.19)THEN
            DONEI = 21
         ELSE IF(SPEC.EQ.73)THEN
            DONEI = 16
         ELSE IF(SPEC.EQ.93)THEN
            DONEI = 41
         ELSE IF(SPEC.EQ.108)THEN
            DONEI = 18
         ELSE IF(SPEC.EQ.119)THEN
            DONEI = 7
         ELSE IF(SPEC.EQ.122)THEN
            DONEi = 32
         ELSE IF(SPEC.EQ.202)THEN
            DONEI = 21
         ELSE IF(SPEC.EQ.242)THEN
            DONEI = 22
         ELSE IF(SPEC.EQ.263 .OR. SPEC.EQ.264)THEN
            DONEI = 14
         ENDIF

     ENDIF

     IF(DONEI.GT.0) THEN
         VOLEQ = EQNUMI(DONEI)
      ELSE IF(DONEF.GT.0) THEN
         VOLEQ = EQNUMF(DONEF)
     ELSE
!           No INGY, find Behre's hyperbola model
        LAST = 53

        FIRST = 1

        DO 7, WHILE (DONEF.EQ.0)
          HALF = (LAST - FIRST +1)/2 + FIRST
          IF(FIA(HALF) .EQ. SPEC)THEN
            DONEF = HALF
          ELSEIF(FIRST .EQ. LAST) THEN
            ERRFLAG = 1
            DONEF = -1
          ELSE IF (FIA(HALF) .LT. SPEC) THEN
            FIRST = HALF
          ELSE
            LAST = HALF - 1
          ENDIF
7 CONTINUE

        IF(DONEF .LT. 0) THEN
          VOLEQ = '616BEHW000'
        ELSE
            IF(FIA(DONEF) .LT. 10) THEN
              WRITE(ASPEC,'(I1)')FIA(DONEF)
              VOLEQ(1:9) = '616BEHW00'
              VOLEQ(10:10) = ASPEC
            ELSE IF (FIA(DONEF) .LT. 100) THEN
              WRITE(ASPEC,'(I2)')FIA(DONEF)
              VOLEQ(1:8) = '616BEHW0'
              VOLEQ(9:10) = ASPEC
            ELSE
              WRITE(ASPEC,'(I3)')FIA(DONEF)
              VOLEQ(1:7) = '616BEHW'
              VOLEQ(8:10) = ASPEC
            ENDIF
        ENDIF
     ENDIF

  ENDIF
  RETURN
  END


!//////////////////////////////////////////////////////////////////
  SUBROUTINE R7_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 VAR,FORST
  INTEGER SPEC,ERRFLAG,FORNUM
  CHARACTER*10 EQNUM(45)
  CHARACTER*10 EQNUMI(60)
  CHARACTER*10 EQNUMC(33),EQNUMF(27),EQNUMD(8)
  INTEGER FIA(41), FIRST, HALF, LAST, DONE,I

!     SPECIES
!     Pacific. silver fir,White fir,       Grand fir,        Subalpine fir,    California red fir,
!     Shasta red fir,     Noble fir,       Port Orford cedar,Alaska cedar,     Western juniper,
!     Subalpine larch,    Western larch,   Incense cedar,    Engelmann spruce, Sitka spruce,
!     Whitebark pine,     Knobcone pine,   Lodgepole pine,   Jeffery pine,     Sugar pine,
!     Western white pine, Ponderosa pine,  Douglas fir,      Redwood,          Pacific yew,
!     Western red cedar,  Hemlock,         Western Hemlock,  Mountain Hemlock, Big leaf maple,
!     Red alder,          White alder,     Pacific madrone,  Golden chinquapin,Oregon ash,
!     Tanoak,             Black cottonwood,Oak species,      Oregon White Oak, Oregon myrtle,
!     Other species,

  DATA (FIA(I),I=1,41)/  11,   15,   17,   19,   20, &
                         21,   22,   41,   42,   64, &
                         72,   73,   81,   93,   98, &
                        101,  103,  108,  116,  117, &
                        119,  122,  202,  211,  231, &
                        242,  260,  263,  264,  312, &
                        351,  352,  361,  431,  542, &
                        631,  747,  800,  815,  981, &
                        999/

  DATA (EQNUM(I),I=1,45)/ &
  'B00BEHW011','B00BEHW015','B00BEHW017','B00BEHW015','B00BEHW021', &
  'B00BEHW021','B00BEHW022','B00BEHW041','B00BEHW042','B00BEHW242', &
  'B00BEHW073','B00BEHW073','B00BEHW081','B00BEHW093','B00BEHW098', &
  'B00BEHW119','B00BEHW108','B00BEHW108','B00BEHW116','B00BEHW117', &
  'B00BEHW119','B00BEHW122','B01BEHW202','B00BEHW211','B00BEHW231', &
  'B00BEHW242','B00BEHW260','B00BEHW260','B00BEHW260','B00BEHW312', &
  'B00BEHW351','B00BEHW361','B00BEHW361','B00BEHW431','B00BEHW542', &
  'B00BEHW631','B00BEHW747','B00BEHW800','B00BEHW800','B00BEHW998', &
  'B00BEHW999', &
  'B02BEHW202','B03BEHW202','B01BEHW202','B00BEHW263'/
!
!  REGION 6 EQUATIONS
!
  DATA (EQNUMI(I),I=1,60)/ &
  'I00FW2W012','I00FW2W017','I00FW2W019','I00FW2W073','I00FW2W093', &
  'I00FW2W108','I00FW2W119','I00FW2W122','I00FW2W202','I00FW2W242', &
  'I00FW2W260','I00FW2W263','I11FW2W012','I11FW2W017','I11FW2W019', &
  'I11FW2W073','I11FW2W093','I11FW2W108','I11FW2W119','I11FW2W122', &
  'I11FW2W202','I11FW2W242','I11FW2W260','I11FW2W263','I12FW2W012', &
  'I12FW2W017','I12FW2W019','I12FW2W073','I12FW2W093','I12FW2W108', &
  'I12FW2W119','I12FW2W122','I12FW2W202','I12FW2W242','I12FW2W260', &
  'I12FW2W263','I13FW2W012','I13FW2W017','I13FW2W019','I13FW2W073', &
  'I13FW2W093','I13FW2W108','I13FW2W119','I13FW2W122','I13FW2W202', &
  'I13FW2W242','I13FW2W260','I13FW2W263','I14FW2W012','I14FW2W017', &
  'I14FW2W019','I14FW2W073','I14FW2W093','I14FW2W108','I14FW2W119', &
  'I14FW2W122','I14FW2W202','I14FW2W242','I14FW2W260','I14FW2W263'/

  DATA (EQNUMF(I),I=1,27)/ &
  'F00FW2W202','F00FW2W242','F00FW2W263','F01FW2W202','F01FW2W242', &
  'F01FW2W263','F02FW2W202','F02FW2W242','F02FW2W263','F03FW2W202', &
  'F03FW2W242','F03FW2W263','F04FW2W202','F04FW2W242','F04FW2W263', &
  'F05FW2W202','F05FW2W242','F05FW2W263','F06FW2W202','F06FW2W242', &
  'F06FW2W263','F07FW2W202','F07FW2W242','F07FW2W263','F08FW2W202', &
  'F08FW2W242','F08FW2W263'/

  DATA (EQNUMC(I),I=1,33)/ &
  'I21FW2W012','I21FW2W017','I21FW2W019','I21FW2W073','I21FW2W093', &
  'I21FW2W108','I21FW2W119','I21FW2W122','I21FW2W202','I21FW2W242', &
  'I21FW2W260','I22FW2W012','I22FW2W017','I22FW2W019','I22FW2W073', &
  'I22FW2W093','I22FW2W108','I22FW2W119','I22FW2W122','I22FW2W202', &
  'I22FW2W242','I22FW2W260','I23FW2W012','I23FW2W017','I23FW2W019', &
  'I23FW2W073','I23FW2W093','I23FW2W108','I23FW2W119','I23FW2W122', &
  'I23FW2W202','I23FW2W242','I23FW2W260'/

  DATA (EQNUMD(I),I=1,8)/ &
  '601DVEW205','601DVEW263','601DVEW015','602DVEW122','602DVEW204', &
  '602DVEW015','602DVEW108','602DVEW122'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    IF(VOLEQ(1:7).EQ.'616BEHW')THEN
!
!  FOUND VALID WESTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF

    IF(VOLEQ(1:7).EQ.'632BEHW')THEN
!
!  FOUND VALID EASTSIDE EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    DO I=1,44
    IF(VOLEQ.EQ.EQNUM(I))THEN
!
!  FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
    DO I=1,60
    IF(VOLEQ.EQ.EQNUMI(I))THEN
!
!  FOUND VALID INGY EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,27
    IF(VOLEQ.EQ.EQNUMF(I))THEN
!
!  FOUND VALID WESTSIDE FLEWELLING EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,33
    IF(VOLEQ.EQ.EQNUMC(I))THEN
!
!  FOUND VALID INGY EQUATION NUMBER - CANADIAN MODEL
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO

    DO I=1,8
    IF(VOLEQ.EQ.EQNUMD(I))THEN
!
!  FOUND VALID DIRECT VOLUME ESTIMATORS
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  DONE = 0

  LAST = 41

  FIRST = 1
  IF(SPEC.EQ.202) THEN
      IF(VAR.EQ.'WC' .OR. VAR.EQ.'wc')THEN
        DONE=44
      ELSEIF((VAR.EQ.'PN') .OR. (VAR.EQ.'pn') .OR. &
         (VAR.EQ.'NC') .OR. (VAR.EQ.'nc') .OR. &
         (VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
         (VAR.EQ.'OC') .OR. (VAR.EQ.'oc') .OR. &
         (VAR.EQ.'OP') .OR. (VAR.EQ.'op')) THEN
        DONE = 44
        IF(FORNUM.EQ.12)DONE=42
      ELSE
            DONE = 44
      ENDIF
  ELSEIF(SPEC.EQ.41)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=13
     ENDIF
  ELSEIF(SPEC.EQ.263)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'PN') .OR. (VAR.EQ.'pn') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc') .OR. &
        (VAR.EQ.'OP') .OR. (VAR.EQ.'op')) THEN
        DONE=45
     ENDIF
  ELSEIF((SPEC.EQ.109).OR.(SPEC.EQ.113).OR.(SPEC.EQ.124).OR. &
     (SPEC.EQ.127))THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=18
     ENDIF
  ELSEIF(SPEC.EQ.92)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=14
     ENDIF
  ELSEIF(SPEC.EQ.212)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=24
     ENDIF
  ELSEIF((SPEC.EQ.801).OR.(SPEC.EQ.805).OR.(SPEC.EQ.807).OR. &
         (SPEC.EQ.811).OR.(SPEC.EQ.818).OR.(SPEC.EQ.821).OR. &
         (SPEC.EQ.839).OR.(SPEC.EQ.333).OR.(SPEC.EQ.730))THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=38
     ELSEIF((VAR.EQ.'NC').OR.(VAR.EQ.'nc'))THEN
        IF(SPEC.EQ.818)DONE=38
     ENDIF
  ELSEIF(SPEC.EQ.542)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=30
     ENDIF
  ELSEIF(SPEC.EQ.251)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
       DONE=25
     ENDIF
  ELSEIF(SPEC.EQ.981)THEN
     IF((VAR.EQ.'CA') .OR. (VAR.EQ.'ca') .OR. &
        (VAR.EQ.'OC') .OR. (VAR.EQ.'oc'))THEN
        DONE=36
     ENDIF
   ENDIF
!
!  FINDS SPECIES INDEX TO FIA ARRAY
!
  DO 5, WHILE (DONE.EQ.0)
     HALF = (LAST - FIRST +1)/2 + FIRST
      IF(FIA(HALF) .EQ. SPEC)THEN
         DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
         ERRFLAG = 1
         DONE = -1
      ELSE IF (FIA(HALF) .LT. SPEC) THEN
         FIRST = HALF
      ELSE
         LAST = HALF - 1
      ENDIF
5 CONTINUE
!
  IF(DONE .LT. 0) THEN
      VOLEQ = 'B00BEHW999'
  ELSE
      VOLEQ=EQNUM(DONE)
  ENDIF
  RETURN
  END
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R8_BEQN(FORST,DIST,SPEC,PROD,VAR,VOLEQ,ERRFLAG)
  CHARACTER*2 GEOAREA,GEOCODES(33)
  CHARACTER*2 PROD,VAR,FORST,DIST
  CHARACTER*10 VOLEQ,VEQTEM
  INTEGER SPEC,ERRFLAG,FORNUM,DISTNUM,FIRST,HALF,LAST,DONE,I,J
  CHARACTER*3 SNSP(92)
  INTEGER SNFIA(92)

!     match species to valid species equation code
  DATA (SNFIA(I),I=1,92)/ &
    10,  57,  90, 107, 110, 111, 115, 121, 123, 126, 128, &
   129, 130, 131, 132, 221, 222, 260, 299, 311, 313, &
   316, 317, 318, 330, 370, 372, 391, 400, 450, 460, &
   471, 491, 521, 531, 540, 541, 543, 544, 552, 555, &
   580, 591, 601, 602, 611, 621, 650, 651, 652, 653, &
   654, 660, 680, 691, 693, 694, 701, 711, 721, 731, &
   740, 743, 762, 802, 806, 812, 813, 819, 820, 822, &
   824, 825, 826, 827, 830, 832, 833, 834, 835, 837, 838, &
   901, 920, 931, 950, 970, 971, 972, 975, 998, 999/

  DATA (SNSP(I), I=1,92)/ &
  '261','100','115','132','110','111','115','121','126','126', &
  '128', &
  '129','132','131','132','221','222','261','132','500','500', &
  '316','300','500','330','370','370','370','400','300','460', &
  '300','300','500','531','541','541','300','544','500','300', &
  '300','300','500','500','611','621','652','300','652','653', &
  '300','300','300','300','693','694','500','300','300','731', &
  '300','300','300','802','806','812','813','800','800','822', &
  '800','800','800','827','827','832','833','800','835','800', &
  '835', &
  '901','300','300','300','970','970','970','970','300','300'/

  DATA (GEOCODES(I), I=1,33)/ &
    '01','02','03','04','05','06','07','08','09','10', &
    '11','12','13','14','15','16','17','18','19','20', &
    '21','22','23','24','25','26','27','28','29','30', &
    '31','32','33'/
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    VEQTEM(1:1)='8'
    VEQTEM(4:7)='DVEE'
    DO I=1,33
    VEQTEM(2:3)=GEOCODES(I)
!  SN
    DO J=1,92
    VEQTEM(8:10)=SNSP(J)
    IF(VOLEQ.EQ.VEQTEM)THEN
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
    ENDDO
  RETURN
  ENDIF

  READ(FORST,'(I2)')FORNUM
  READ(DIST,'(I2)')DISTNUM

  IF (FORNUM.EQ.1)THEN
!              // alabama
          IF (DISTNUM.EQ.1) THEN
                 GEOAREA = '13';
          ELSE IF (DISTNUM.EQ.3) THEN
                 GEOAREA = '15';
          ELSE IF (DISTNUM.EQ.4) THEN
                 GEOAREA = '16';
          ELSE IF (DISTNUM.EQ.5 .OR. DISTNUM.EQ.6) THEN
             GEOAREA = '14';
          ELSE
             GEOAREA = '17';
     ENDIF

   ELSE IF (FORNUM.EQ.2)THEN
!                     // daniel boone
                 GEOAREA = '10';

   ELSE IF (FORNUM.EQ.3)THEN
!                     // Chattahoochee/Oconee
                 GEOAREA = '01';
                 IF (DISTNUM .EQ. 8) GEOAREA = '06';

   ELSE IF (FORNUM.EQ.4 .OR. FORNUM .EQ. 60)THEN
!                     // Cherokee
!                      // land between the lakes
                 GEOAREA = '01';

   ELSE IF (FORNUM.EQ.5)THEN
!                     // Florida
                 GEOAREA = '03';

   ELSE IF (FORNUM.EQ.6) THEN
!                     // Kisatchie
                 GEOAREA = '02';
                 IF (DISTNUM .EQ. 6) GEOAREA = '09';

   ELSE IF (FORNUM .EQ. 7) THEN
!                     // Mississippi

          IF (DISTNUM.EQ.1) THEN
                 GEOAREA = '19';
          ELSE IF (DISTNUM.EQ.2) THEN
                 GEOAREA = '20';
          ELSE IF (DISTNUM.EQ.4) THEN
                 GEOAREA = '21';
          ELSE IF (DISTNUM.EQ.5) THEN
             GEOAREA = '22';
          ELSE IF (DISTNUM.EQ.6) THEN
             GEOAREA = '18';
          ELSE
             GEOAREA = '23';
     END IF

   ELSE IF (FORNUM .EQ. 8) THEN
!                     // GW/Jeff
      GEOAREA = '11';

     IF (DISTNUM.EQ.11 .OR. DISTNUM.EQ.12 .OR. DISTNUM.EQ.13 .OR. &
         DISTNUM.EQ.14 .OR. DISTNUM.EQ.15 .OR. DISTNUM.EQ.16) THEN
        GEOAREA = '12';
     END IF

   ELSE IF (FORNUM .EQ. 9) THEN
!                     // Ouachita
          GEOAREA = '31';
          IF (DISTNUM.EQ.1 .OR. DISTNUM.EQ.6) THEN
                 GEOAREA = '30';
          ELSE IF (DISTNUM.EQ.12) THEN
             GEOAREA = '32';
     END IF

  ELSE IF (FORNUM .EQ. 10) THEN
!                     // Ozark/St Francis
                 GEOAREA = '04';
             IF (DISTNUM .EQ. 7) GEOAREA = '05';

   ELSE IF (FORNUM .EQ. 11) THEN
!                     // North Carolina
                 GEOAREA = '01';
                 IF (DISTNUM .EQ. 3)THEN
                    GEOAREA = '07';
                 ELSE IF (DISTNUM .EQ. 10) THEN
                        GEOAREA = '08';
                 ENDIF

   ELSE IF (FORNUM .EQ. 12)THEN
!                     // Francis Marion/Sumpter
             GEOAREA = '24';
                 IF(DISTNUM .EQ. 2) THEN
                    GEOAREA = '01';
                 ELSE IF(DISTNUM .EQ. 5)THEN
                        GEOAREA = '25';
                 ENDIF

   ELSE IF (FORNUM.EQ.13) THEN
!                      // Texas
                 IF(DISTNUM .EQ. 1) THEN
                    GEOAREA = '26';
                 ELSE IF(DISTNUM .EQ. 3)THEN
                        GEOAREA = '27';
                 ELSE IF(DISTNUM .EQ. 4)THEN
                        GEOAREA = '29';
                 ELSE
                        GEOAREA = '28';
                 ENDIF

   ELSE IF (FORNUM.EQ.36)THEN
!            // savannah river
                 GEOAREA = '25';

  ENDIF

!     CREATE THE VOLUME EQUATION NUMBER

  VOLEQ(1:1) = '8'
  VOLEQ(2:3) = GEOAREA
  VOLEQ(4:7) = 'DVEE'

!     FIND CORRECT SPECIES
  DONE = 0
  FIRST = 1
  LAST = 92
  DO 5, WHILE (DONE.EQ.0)
     HALF = (LAST - FIRST +1)/2 + FIRST
      IF(SNFIA(HALF) .EQ. SPEC)THEN
         DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
         ERRFLAG = 1
         DONE = -1
     ELSE IF (SNFIA(HALF) .LT. SPEC) THEN
         FIRST = HALF
      ELSE
         LAST = HALF - 1
      ENDIF
5 CONTINUE
  IF(DONE .LT. 0) THEN
    IF(SPEC.LT.300) THEN
! USE OTHER SOFTWOOD 299
      DONE = 19
    ELSE
      DONE = 92
    ENDIF
  ENDIF

  VOLEQ(8:10) = SNSP(DONE)
  RETURN
  END
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R8_CEQN(FORST,DIST,SPEC,PROD,VAR,VOLEQ,ERRFLAG)
  CHARACTER*1 GEOAREA,TOPCODE(5),ICHAR
  CHARACTER*2 PROD,VAR,FORST,DIST
  CHARACTER*10 VOLEQ,VEQTEM
  CHARACTER*3 SNSP(110)
  INTEGER SNFIA(110)
  INTEGER SPEC,ERRFLAG,FORNUM,DISTNUM,FIRST,HALF,LAST,DONE,I,J,K

!     match species to valid species equation code
!Updated the species list to include R8 cruise species (2025/05/23)
!      DATA (SNFIA(I),I=1,92)/
!     >  10,  57,  90, 107, 110, 111, 115, 121, 123, 126, 128,
!     > 129, 130, 131, 132, 221, 222, 260, 299, 311, 313,
!     > 316, 317, 318, 330, 370, 372, 391, 400, 450, 460,
!     > 471, 491, 521, 531, 540, 541, 543, 544, 552, 555,
!     > 580, 591, 601, 602, 611, 621, 650, 651, 652, 653,
!     > 654, 660, 680, 691, 693, 694, 701, 711, 721, 731,
!     > 740, 743, 762, 802, 806, 812, 813, 819, 820, 822,
!     > 824, 825, 826, 827, 830, 832, 833, 834, 835, 837, 838,
!     > 901, 920, 931, 950, 970, 971, 972, 975, 998, 999/

!      DATA (SNSP(I), I=1,92)/
!     &'261','100','115','132','110','111','115','121','126','126',
!     &'128',
!     &'129','132','131','132','221','222','261','132','500','500',
!     &'316','300','500','330','370','370','370','400','300','460',
!     &'300','300','500','531','541','541','300','544','500','300',
!     &'300','300','500','500','611','621','652','300','652','653',
!     &'300','300','300','300','693','694','500','300','300','731',
!     &'300','300','300','802','806','812','813','800','800','822',
!     &'800','800','800','827','827','832','833','800','835','800',
!     &'835',
!     &'901','300','300','300','970','970','970','970','300','300'/
  DATA (SNFIA(I),I=1,110)/ &
    10, 57, 90,100,107,110,111,115,121,123, &
   126,128,129,130,131,132,197,221,222,260, &
   261,299,300,311,313,314,316,317,318,330, &
   370,372,391,400,404,450,460,471,491,500, &
   521,531,540,541,543,544,545,546,550,552, &
   555,580,591,601,602,611,621,650,651,652, &
   653,654,660,680,691,693,694,701,711,721, &
   731,740,742,743,762,800,802,804,806,812, &
   813,817,819,820,822,823,824,825,826,827, &
   828,830,831,832,833,834,835,837,838,901, &
   920,930,931,950,970,971,972,975,998,999/
  DATA (SNSP(I), I=1,110)/ &
   '261','100','115','100','107','110','111','115','121','123', &
   '126','128','129','132','131','132','197','221','222','261', &
   '261','132','300','500','313','314','316','317','318','330', &
   '370','370','370','400','404','300','460','300','300','500', &
   '521','531','541','541','300','544','545','546','550','500', &
   '300','580','300','601','602','611','621','650','651','652', &
   '653','300','300','300','691','693','694','500','711','300', &
   '731','300','742','300','762','800','802','804','806','812', &
   '813','817','800','820','822','823','800','825','826','827', &
   '828','830','831','832','833','834','835','837','835','901', &
   '920','930','300','950','970','970','970','970','300','300'/
!
  DATA TOPCODE / '1','4','7','8','9' /
!
!  SEARCH FOR VALID EQUATION NUMBER
!
  IF(SPEC.EQ.9999)THEN
    VEQTEM(1:1)='8'
    VEQTEM(4:7)='CLKE'
    DO I=1,7
    WRITE(ICHAR,'(I1)')I
    VEQTEM(2:2)=ICHAR
    DO J=1,5
    VEQTEM(3:3)=TOPCODE(J)
!  SN
    DO K=1,110
    VEQTEM(8:10)=SNSP(K)
    IF(VOLEQ.EQ.VEQTEM)THEN
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
    ENDDO
    ENDDO
  RETURN
  ENDIF

!     FIND CORRECT FOREST
  READ(FORST,'(I2)')FORNUM
  READ(DIST,'(I2)')DISTNUM

!     SET GEOCODE
!     SET DIAMETER BY PRODUCT

  IF (FORNUM.EQ.1)THEN
!             // alabama
          GEOAREA = '4';
          IF (DISTNUM.EQ.3) GEOAREA = '1'

   ELSE IF (FORNUM.EQ.2 .OR. FORNUM.EQ.4 .OR. FORNUM .EQ. 8 .OR. &
            FORNUM .EQ. 60)THEN
!                     // daniel boone
!                     // Cherokee
!                     // GW/Jeff
!                 // land between the lakes
                 GEOAREA = '3';

   ELSE IF (FORNUM.EQ.3)THEN
!                     // Chattahoochee/Oconee
                 GEOAREA = '3';
                 IF (DISTNUM .EQ. 8) GEOAREA = '2';

   ELSE IF (FORNUM.EQ.5 .OR. FORNUM.EQ.36)THEN
!                     // Florida
!           // savannah river
                 GEOAREA = '1';

   ELSE IF (FORNUM.EQ.6 .OR. FORNUM.EQ.13) THEN
!                     // Kisatchie
!                     // Texas
                 GEOAREA = '5';

   ELSE IF (FORNUM .EQ. 7) THEN
!                     // Mississippi
                 GEOAREA = '5';

        IF (DISTNUM .EQ. 6) THEN
                 GEOAREA = '7';
        ELSE IF (DISTNUM.EQ.7 .OR. DISTNUM.EQ.17) THEN
                 GEOAREA = '4';
        ENDIF
   ELSE IF (FORNUM .EQ. 9) THEN
!                     // Ouachita
                 GEOAREA = '6';
   ELSE IF (FORNUM .EQ. 10) THEN
!                     // Ozark/St Francis
                 GEOAREA = '6';
             IF (DISTNUM .EQ. 7) GEOAREA = '7';

   ELSE IF (FORNUM .EQ. 11) THEN
!              // North Carolina
           GEOAREA = '3';
          IF (DISTNUM .EQ. 3)THEN
               GEOAREA = '1';
          ELSE IF (DISTNUM .EQ. 10) THEN
               GEOAREA = '2';
          ENDIF

   ELSE IF (FORNUM .EQ. 12)THEN
!                     // Francis Marion/Sumpter
             GEOAREA = '2';
             IF(DISTNUM .EQ. 2) THEN
                GEOAREA = '3';
             ELSE IF(DISTNUM .EQ. 5)THEN
                GEOAREA = '1';
             ENDIF
  ENDIF

!     CREATE THE VOLUME EQUATION NUMBER
  VOLEQ(1:1) = '8'
  VOLEQ(2:2) = GEOAREA
  !Comment out for PROD = '08' (2025/05/23)
!       IF (PROD.EQ.'08') THEN
!C           USE PRODUCT 08 LOGIC
!          VOLEQ(3:3) = '8'
!       ELSE
!         VOLEQ(3:3) = '1'
!       ENDIF
  VOLEQ(3:3) = '1'
  VOLEQ(4:7) = 'CLKE'

!     FIND CORRECT SPECIES
  DONE = 0
  FIRST = 1
  LAST = 110
  DO 5, WHILE (DONE.EQ.0)
     HALF = (LAST - FIRST +1)/2 + FIRST
      IF(SNFIA(HALF) .EQ. SPEC)THEN
         DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
         ERRFLAG = 1
         DONE = -1
     ELSE IF (SNFIA(HALF) .LT. SPEC) THEN
         FIRST = HALF
      ELSE
         LAST = HALF - 1
      ENDIF
5 CONTINUE
  IF(DONE .LT. 0) THEN
    IF(SPEC.LT.300) THEN
! USE OTHER SOFTWOOD 299
      DONE = 22
    ELSE
      DONE = 110
    ENDIF
  ENDIF

  VOLEQ(8:10) = SNSP(DONE)

   RETURN
  END


!//////////////////////////////////////////////////////////////////
  SUBROUTINE R9_EQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ,VEQTEM
  CHARACTER*2 GEOAREA,VAR, FORST,GEOCODES(6)
  CHARACTER*3 LSSP(69),CSSP(97),NESP(108),SNSP(92),ASPEC
  INTEGER SPEC,ERRFLAG,FORNUM,LSFIA(69),CSFIA(97),NEFIA(108), &
           SNFIA(92)
  INTEGER FIRST, HALF, LAST, DONE, I, J

  DATA (LSSP(I), I=1,69)/ &
  '012','068','071','091','094','095','105','125','125','129', &
  '130','241','261','001','313','314','315','316','317','318', &
  '319','371','375','391','402','403','407','421','462','491', &
  '500','531','541','543','544','601','602','660','693','701', &
  '731','741','742','743','746','766','761','762','763','766', &
  '802','804','809','823','826','833','837','901','920','922', &
  '923','931','935','951','972','975','977','300','999'/

  DATA (CSSP(I), I=1,97)/ &
  '068','068','110','129','131','132','221','001','316','316', &
  '316','318','331','373','391','400','401','402','403','404', &
  '405','407','408','409','400','450','460','460','471','490', &
  '500','521','531','541','541','543','541','543','543','552', &
  '571','601','602','611','621','641','651','653','680','690', &
  '691','693','694','701','711','731','742','741','742','746', &
  '746','762','802','823','806','812','813','813','813','823', &
  '824','823','823','816','813','830','813','832','833','813', &
  '835','835','837','901','920','920','931','951','970','970', &
  '970','970','970','970','999','999','999'/


  DATA (NESP(I), I=1,108)/ &
  '012','068','068','068','071','097','097','094','095','097', &
  '105','105','110','105','125','105','105','129','105','131', &
  '125','068','261','261','105','313','318','315','316','318', &
  '318','330','332','341','355','371','371','373','374','375', &
  '375','391','400','400','400','400','400','462','490','500', &
  '521','531','541','541','543','541','541','591','601','602', &
  '611','621','641','650','651','653','660','691','693','701', &
  '711','712','731','741','742','746','742','746','760','761', &
  '762','802','802','832','806','833','813','806','823','832', &
  '802','806','830','831','832','833','835','837','901','920', &
  '931','951','951','970','970','970','999','004'/

  DATA (SNSP(I), I=1,92)/ &
  '261','100','115','132','110','111','115','121','126','126', &
  '128', &
  '129','132','131','132','221','222','261','132','500','500', &
  '316','300','500','330','370','370','370','400','300','460', &
  '300','300','500','531','541','541','300','544','500','300', &
  '300','300','500','500','611','621','652','300','652','653', &
  '300','300','300','300','693','694','500','300','300','731', &
  '300','300','300','802','806','812','813','800','800','822', &
  '800','800','800','827','827','832','833','800','835','800', &
  '835', &
  '901','300','300','300','970','970','970','970','300','300'/

  DATA (LSFIA(I), I=1,69)/ &
   12,  68,  71,   91,  94,  95, 105, 125, 125, 129, &
   130, 241, 261, 299, 313, 314, 315, 316, 317, 318, &
   319, 371, 375, 391, 402, 403, 407, 421, 462, 491, &
   500, 531, 541, 543, 544, 601, 602, 660, 693, 701, &
   731, 741, 742, 743, 746, 760, 761, 762, 763, 766, &
   802, 804, 809, 823, 826, 833, 837, 901, 920, 922, &
   923, 931, 935, 951, 972, 975, 977, 993, 994/

  DATA (CSFIA(I), I =1,97)/ &
    57,  68, 110, 129, 131, 132, 221, 299, 313, 316, &
   317, 318, 331, 373, 391, 400, 401, 402, 403, 404, &
   405, 407, 408, 409, 410, 450, 461, 462, 471, 491, &
   500, 521, 531, 540, 541, 543, 544, 545, 546, 552, &
   571, 601, 602, 611, 621, 641, 651, 653, 680, 690, &
   691, 693, 694, 701, 711, 731, 740, 741, 742, 743, &
   746, 762, 802, 804, 806, 812, 813, 817, 822, 823, &
   824, 825, 826, 827, 828, 830, 831, 832, 833, 834, &
   835, 836, 837, 901, 920, 922, 931, 951, 970, 971, &
   972, 974, 975, 977, 991, 992, 994/

  DATA (NEFIA(I), I=1,108)/ &
    12,  43,  57,  68,  71,  90,  91,  94,  95,  97, &
   100, 105, 110, 123, 125, 126, 128, 129, 130, 131, &
   132, 241, 260, 261, 299, 313, 314, 315, 316, 317, &
   318, 330, 332, 341, 356, 371, 372, 373, 374, 375, &
   379, 391, 400, 403, 405, 407, 409, 462, 491, 500, &
   521, 531, 540, 541, 543, 544, 545, 591, 601, 602, &
   611, 621, 641, 650, 651, 653, 660, 691, 693, 701, &
   711, 712, 731, 741, 742, 743, 744, 746, 760, 761, &
   762, 800, 802, 804, 806, 812, 813, 817, 823, 825, &
   826, 827, 830, 831, 832, 833, 835, 837, 901, 922, &
   931, 951, 952, 970, 972, 975, 994, 998/

  DATA (SNFIA(I),I=1,92)/ &
    10,  57,  90, 107, 110, 111, 115, 121, 123, 126, 128, &
   129, 130, 131, 132, 221, 222, 260, 299, 311, 313, &
   316, 317, 318, 330, 370, 372, 391, 400, 450, 460, &
   471, 491, 521, 531, 540, 541, 543, 544, 552, 555, &
   580, 591, 601, 602, 611, 621, 650, 651, 652, 653, &
   654, 660, 680, 691, 693, 694, 701, 711, 721, 731, &
   740, 743, 762, 802, 806, 812, 813, 819, 820, 822, &
   824, 825, 826, 827, 830, 832, 833, 834, 835, 837, 838, &
   901, 920, 931, 950, 970, 971, 972, 975, 998, 999/
!
!  SEARCH FOR VALID EQUATION NUMBER
!  FIRST, SEARCH FOR CLKE OR DVEE
!
  IF(SPEC.EQ.9999)THEN
    IF((VOLEQ(1:7).EQ.'900CLKE').OR.(VOLEQ(4:7).EQ.'DVEE'))THEN
!  LS
      DO J=1,69
      IF(VOLEQ(8:10).EQ.LSSP(J))THEN
        SPEC=8888
        RETURN
      ENDIF
      ENDDO
!  CS
      DO J=1,97
      IF(VOLEQ(8:10).EQ.CSSP(J))THEN
        SPEC=8888
        RETURN
      ENDIF
      ENDDO
!  NE
      DO J=1,108
      IF(VOLEQ(8:10).EQ.NESP(J))THEN
        SPEC=8888
        RETURN
      ENDIF
      ENDDO
!  SN
      DO J=1,92
      IF(VOLEQ(8:10).EQ.SNSP(J))THEN
        SPEC=8888
        RETURN
      ENDIF
      ENDDO
      RETURN
    ELSE
!     NOT A VALID REGION 9 EQUATION
      RETURN
    ENDIF
  ENDIF
!
  IF(VOLEQ(1:7).EQ.'900CLKE')THEN
!     NEW CLARK'S PROFILE MODEL VOLUME EQUATION NUMBERS
!     MAKE SURE SPEC IS A 3 CHARACTER FIELD.
!
    WRITE(ASPEC,'(I3)')SPEC

    IF(SPEC .LT.10)THEN
      ASPEC(1:2) = '00'
    ELSEIF(SPEC.LT.100)THEN
      ASPEC(1:1) = '0'
    ENDIF

    VOLEQ(8:10) = ASPEC

    RETURN
  ELSEIF(VOLEQ(4:7).EQ.'DVEE')THEN
!
!     DIRECT VOLUME ESTIMATORS
!     FIND CORRECT SPECIES
    FIRST = 1
    DONE = 0
    IF(VAR.EQ.'LS' .OR. VAR.EQ.'ls') THEN
      LAST = 69
      DO 5, WHILE (DONE.EQ.0)
      HALF = (LAST - FIRST +1)/2 + FIRST
      IF(LSFIA(HALF) .EQ. SPEC)THEN
        DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
        ERRFLAG = 1
        DONE = -1
      ELSEIF(LSFIA(HALF) .LT. SPEC) THEN
        FIRST = HALF
      ELSE
        LAST = HALF - 1
      ENDIF
5 CONTINUE
      IF(DONE .LT. 0) DONE = 69

      VOLEQ(8:10) = LSSP(DONE)
    ELSE IF(VAR.EQ.'CS' .OR. VAR.EQ.'cs')THEN
      LAST = 97
      DO 15, WHILE (DONE.EQ.0)
      HALF = (LAST - FIRST +1)/2 + FIRST
      IF(CSFIA(HALF) .EQ. SPEC)THEN
        DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
        ERRFLAG = 1
        DONE = -1
      ELSE IF (CSFIA(HALF) .LT. SPEC) THEN
        FIRST = HALF
      ELSE
        LAST = HALF - 1
      ENDIF
15 CONTINUE
      IF(DONE .LT. 0) DONE = 97
       VOLEQ(8:10) = CSSP(DONE)
    ELSE IF(VAR.EQ.'NE' .OR. VAR.EQ.'ne')THEN
      LAST = 108
      DO 25, WHILE (DONE.EQ.0)
      HALF = (LAST - FIRST +1)/2 + FIRST
      IF(NEFIA(HALF) .EQ. SPEC)THEN
        DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
        ERRFLAG = 1
        DONE = -1
      ELSE IF (NEFIA(HALF) .LT. SPEC) THEN
        FIRST = HALF
      ELSE
        LAST = HALF - 1
      ENDIF
25 CONTINUE
      IF(DONE .LT. 0) DONE = 108
      VOLEQ(8:10) = NESP(DONE)
    ELSE
!
!  SOUTHERN VARIANT RUNNING REGION 9 FORESTS AND GEVORKIANTZ METH=5
!
      LAST = 92
      DO 35, WHILE (DONE.EQ.0)
      HALF = (LAST - FIRST +1)/2 + FIRST
      IF(SNFIA(HALF) .EQ. SPEC)THEN
        DONE = HALF
      ELSEIF(FIRST .EQ. LAST) THEN
        ERRFLAG = 1
        DONE = -1
      ELSE IF (SNFIA(HALF) .LT. SPEC) THEN
        FIRST = HALF
      ELSE
        LAST = HALF - 1
      ENDIF
35 CONTINUE
      IF(DONE .LT. 0) DONE = 92
      VOLEQ(8:10) = SNSP(DONE)
      ENDIF
    RETURN
  ENDIF
!     End FVS check equation
!     NEW CLARK'S PROFILE MODEL VOLUME EQUATION NUMBERS
!     MAKE SURE SPEC IS A 3 CHARACTER FIELD.
!
    VOLEQ(1:7)='900CLKE'
    WRITE(ASPEC,'(I3)')SPEC

    IF(SPEC .LT.10)THEN
      ASPEC(1:2) = '00'
    ELSEIF(SPEC.LT.100)THEN
      ASPEC(1:1) = '0'
    ENDIF

    VOLEQ(8:10) = ASPEC

  RETURN
  END
!//////////////////////////////////////////////////////////////////
  SUBROUTINE R10_EQN(FORST,SPEC,VOLEQ,ERRFLAG)
  CHARACTER*10 VOLEQ
  CHARACTER*2 FORST
  INTEGER SPEC,ERRFLAG,FORNUM
  CHARACTER*10 TONEQN(27),CHUEQN(27),OTHEREQN(27)
  INTEGER FIA(27), FIRST, HALF, LAST, DONE,I
!
!     AK SPECIES LIST
!     Pacific silver fir  Subalpine fir     Alaska yellow Cedar  Tamarack
!     White spruce        Lutz's spruce     Black spruce         Sitka spurce
!     Lodgepole pine      Western redcedar  Western hemlock      Mountain Hemlock
!     Other Softwood      Alder species     Red alder            Paper birch
!     Alaska birch        Balsam poplar     Quaking aspen        Black cottonwood
!     Willow species      Scouler's willow  Other hardwood
!
!     LUTZ SPRUCE IS TREATED AS FIA CODE 94 (WHITE SPRUCE)
!     Added second growth equation for 042,098, 242,263 using species code 2042,2098,2242 and 2263 (2024/10/07)
  DATA FIA / 011, 019, 042, 071, 094, &
             094, 095, 098, 108, 242, &
             263, 264, 299, 350, 351, &
             375, 376, 741, 746, 747, &
             920, 928, 998,2042,2098, &
            2242,2263 /
!
!     DEFAULT VOLUME EQUATIONS FOR LOCATION CODES: 703, 1002, 1003, 1005,
!     8134, 8135, 8112
  DATA TONEQN / &
   'A00F32W260','A00F32W260','A00F32W042','A00DVEW094','A00DVEW094', &
   'A00DVEW094','A00DVEW094','A00F32W098','A00F32W260','A00F32W242', &
   'A00F32W260','A00F32W260','A00DVEW094','A32CURW351','A32CURW351', &
   'A00DVEW375','A00DVEW375','A00DVEW747','A00DVEW375','A00DVEW747', &
   'A00DVEW747','A00DVEW747','A00DVEW747','A02F32W042','A02F32W098', &
   'A02F32W242','A02F32W263'/

!     DEFAULT VOLUME EQUATIONS FOR LOCATION CODES: 1004, 713, 720, 7400, 7401,
!     7402, 7403, 7404, 7405, 7406, 7407, 7408
  DATA CHUEQN / &
   'A01DEMW000','A01DEMW000','A00DVEW094','A00DVEW094','A00DVEW094', &
   'A00DVEW094','A00DVEW094','A00F32W098','A01DEMW000','A01DEMW000', &
   'A00F32W260','A01DEMW000','A00DVEW094','A32CURW351','A32CURW351', &
   'A00DVEW375','A00DVEW375','A00DVEW747','A00DVEW375','A00DVEW747', &
   'A00DVEW747','A00DVEW747','A00DVEW747','A02F32W042','A02F32W098', &
   'A02F32W242','A02F32W263'/

!     OTHER VALID EQUATIONS THAT CAN BE USED -- THESE ARE NOT DEFAULT
!     EQUATION NUMBERS.
  DATA OTHEREQN / &
   'A00FW2W042','A16DEMW042','A61DEMW042','A32DEMW042','A01DVEW094', &
   'A02DVEW094','A16DEMW098','A00FW2W098','A02F32W098','A02FW2W098', &
   'A02DEMW000','A32CURW000','A32DEMW098','A61DEMW098','A00FW2W242', &
   'A61DEMW242','A16DEMW242','A32DEMW242','A00FW2W260','A02F32W260', &
   'A02FW2W260','A16CURW260','A01DVEW375','A01DVEW747','A00DVEW108', &
   'A00DVEW310','A00DVEW351'/
!
!     SEARCH FOR VALID EQUATION NUMBER IN TONEQN AND CHUEQN
  IF(SPEC.EQ.9999)THEN
    DO I=1,27
    IF((VOLEQ.EQ.TONEQN(I)).OR.(VOLEQ.EQ.CHUEQN(I)))THEN
!
!     FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
!
!     SEARCH FOR VALID EQUATION NUMBER IN OTHEREQN
    DO I=1,27
    IF(VOLEQ.EQ.OTHEREQN(I))THEN
!
!     FOUND VALID EQUATION NUMBER
!
      SPEC=8888
      RETURN
    ENDIF
    ENDDO
  RETURN
  ENDIF
!
  READ(FORST,'(I2)')FORNUM
  LAST = 27
  FIRST = 1
  DONE = 0
!
!     CHECK FOR VALID SPECIES
  DO WHILE (DONE.EQ.0)
     HALF = (LAST - FIRST +1)/2 + FIRST
      IF(FIA(HALF) .EQ. SPEC)THEN
         DONE = HALF
      ELSE IF(FIRST .EQ. LAST) THEN
         ERRFLAG = 1
         DONE = -1
      ELSE IF (FIA(HALF) .LT. SPEC) THEN
         FIRST = HALF
      ELSE
         LAST = HALF - 1
      ENDIF
  ENDDO

!     IF SPECIES NOT FOUND, USE OTHER SOFTWOOD
!     LUTZ SPRUCE WILL NEVER BE FOUND IN THE ABOVE LOOP SINCE
!     IT DOES NOT HAVE A FIA CODE. AS SUCH IT WILL BE ASSIGNED
!     THE EQUATION NUMBER FOR OTHER SOFTWOOD (INDEX 13).
  IF(DONE .LT. 0) DONE = 13
!
!     DETERMINE DEFAULT EQUATION NUMBER TO USE BASED ON INDEX (DONE)
!     AND FOREST NUMBER
  IF(FORNUM .EQ. 4) THEN
     VOLEQ = CHUEQN(DONE)
  ELSE
     VOLEQ = TONEQN(DONE)
  ENDIF

  RETURN
  END

! =======================================================================
  SUBROUTINE FIAVOLEQDEF(VAR,REGN,FORST,DIST,SPEC,VOLEQ,ERRFLAG)

!    SUBROUTINE WILL RETURN THE DEFAULT FIA VOLUME EQUATION NUMBER
!        SPEC = 3 DIGIT FIA SPECIES CODE
! ----------------------------------------------------------------------
! FIA species default equation
! species     common name           WOR     WWA    EOR    EWA     CA
! ----------  -------------------   ---     ---    ---    ---     --
! 11          Pacific silver fir     11     11     10     10     --
! 14          Bristlecone fir        --     --     --     --     18
! 15          White fir              23     --     10     --     23
! 17          Grand fir              11     11     10     10     23
! 19          Subalpine fir          11     11     10     10     18
! 20          California red fir     --     --     --     --     18
! 21          Shasta red fir         18     18     18     18     18
! 22          Noble fir              11     11     10     10     18
! 41          Port-Orford-cedar      19     19     19     19      8
! 42          Alaska-cedar            9      9      8      8      8
! 50          Cypress                --     --     --     --     19
! 51          Arizona cypress        --     --     --     --     19
! 56          McNabb cypress         --     --     --     --     19
! 62          California juniper     --     --     --     --     21
! 64          Western juniper        21     21     21     21     21
! 65          Utah juniper           --     --     --     --     21
! 72          Subalpine larch        --     22     --     22     --
! 73          Western larch          --     22     22     22     --
! 81          Incense cedar          19     19     19     19     19
! 92          Brewer spruce          13     --     13     --     12
! 93          Engelmann spruce       13     13     12     12     12
! 98          Sitka spruce           13     13     --     --     12
! 101          Whitebark pine        15     15     15     15     20
! 102          Bristlecone pine      --     --     --     --     16
! 103          Knobcone pine         15     --     15     --     16
! 104          Foxtail pine          --     --     --     --     16
! 108          Lodgepole pine        15     15     15     15     16
! 109          Coulter pine          --     --     --     --     5
! 113          Limber pine           15     --     15     --     16
! 116          Jeffrey pine           5     --      4     --     5
! 117          Sugar pine            20     20     20     20     20
! 119          Western white pine    15     15     15     15     20
! 120          Bishop pine           --     --     --     --     16
! 122          Ponderosa pine(>=5dbh) 5      4      4      4      5
! 122          Ponderosa pine(<5�dbh) 5      5      5      5      5
! 124          Monterey pine         --     --     --     --     16
! 127          Gray pine             --     --     --     --     5
! 130          Scotch pine           17     17     17     17     17
! 133          Singleleaf pinyon     --     --     --     --     21
! 137          Washoe pine           --     --     --     --     5
! 201          Bigcone Douglas-fir   --     --     --     --     3
! 202          Douglas-fir            1      1      2      2     3
! 211          Redwood               24     --     --     --     24
! 212          Giant Sequoia         24     --     --     --     24
! 231          Pacific yew            9      9      8      8     8
! 242          Western redcedar       9      9      8      8     8
! 251          California nutmeg     --     --     --     --     8
! 263          Western hemlock        6      6      6      6     6
! 264          Mountain hemlock      17     17     17     17     17
! 299          Unknown Conifer       17     17     17     17     17
! 312          Bigleaf maple         37     26     37     26     37
! 313          Boxelder              --     --     --     --     38
! 321          Rocky Mountain maple  --     --     --     --     --
! 322          Bigtooth maple        --     --     --     --     --
! 330          California buckeye    --     --     --     --     43
! 341          Tree of heaven        --     --     --     --     26
! 351          Red alder             26     25     26     25     26
! 352          White alder           26     --     26     --     26
! 361          Pacific madrone       40     26     40     26     40
! 374          Water birch           --     --     --     --     26
! 375          Paper birch           --     --     --     --     --
! 376          Western paper birch   --     26     --     26     --
! 431          Golden chinkapin      32     26     --     26     32
! 475 Curlleaf mountain-mahogany     --     --     45     --     45
! 492          Pacific dogwood       --     26     --     26     26
! 500          Hawthorn              --     --     --     --     --
! 510          Eucalyptus            26     --     --     --     31
! 542          Oregon ash            38     26     38     26     38
! 590          Holly                 26     26     26     26     26
! 600          Walnut                26     26     26     --     38
! 631          Tanoak                34     --     --     --     34
! 660          Apple                 26     26     26     26     42
! 730          California sycamore   26     26     26     26     42
! 740          Cottonwood and poplar --     --     --     --     --
! 741          Balsam poplar         --     --     --     --     --
! 742          Eastern cottonwood    --     --     --     --     --
! 745          Plains cottonwood     --     --     --     --     --
! 746          Quaking aspen         26     26     26     26     28
! 747          Black cottonwood      26     26     26     26     27
! 748          Fremont poplar        --     --     --     --     27
! 755          Mesquite              --     --     --     --     --
! 760          Cherry                26     26     26     26     26
! 800          Oak-deciduous         --     --     --     --     43
! 801          California live oak   --     --     --     --     43
! 805          Canyon live oak       42     --     --     --     42
! 807          Blue oak              --     --     --     --     39
! 810          Emory oak             --     --     --     --     --
! 811          Englemann oak         --     --     --     --     36
! 815          Oregon white oak      41     26     41     26     41
! 818          California black oak  38     --     38     26     38
! 821          California white oak  --     --     --     --     35
! 839          Interior live oak     --     --     --     --     44
! 901          Black locust          --     --     --     --     41
! 920          Willow                26     26     26     26     40
! 981          California-laurel     33     --     --     --     33
! 998          Unknown hardwood      26     26     26     26     41
! 999          Unknown Tree          26     26     26     26     41
! ------------------------------------------------------------------------
  CHARACTER*2 FORST,DIST,VAR
  CHARACTER*10 VOLEQ
  INTEGER ERRFLAG,REGN,SPEC

  IF(VAR .EQ. '  ') THEN
     CALL GETVARIANT(REGN,FORST,DIST,VAR)
  ENDIF
  IF(REGN .EQ. 5) THEN
     CALL R5_PNWEQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  ELSE IF(REGN.EQ.6)THEN
     CALL R6_PNWEQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  ELSE
     ERRFLAG=1
  ENDIF

  RETURN
  END

! FIA DEFAULT VOLUME EQUATION FOR CA
  SUBROUTINE R5_PNWEQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  CHARACTER*2 FORST,DIST,VAR
  CHARACTER*10 VOLEQ, EQNUM(91)
  INTEGER ERRFLAG,SPEC,FIA(91)
  INTEGER DONE,HALF,FIRST,LAST

  DATA (FIA(I),I=1,91)/ &
    11, 14, 15, 17, 19, 20, 21, 22, 41, 42, &
    50, 51, 56, 62, 64, 65, 72, 73, 81, 92, &
    93, 98,101,102,103,104,108,109,113,116, &
   117,119,120,122,124,127,130,133,137,142, &
   201,202,211,212,231,242,251,263,264,299, &
   301,312,313,330,333,341,351,352,361,374, &
   431,475,478,492,510,542,590,600,631,660, &
   671,730,746,747,748,760,768,800,801,805, &
   807,811,815,818,821,839,901,920,981,998, &
   999/
  DATA (EQNUM(I),I=1,91)/ &
!        11,14, 15, 17, 19,
  '616TRFW264','516TRFW021','516TRFW015','516TRFW015','516TRFW021', &
!        20,21, 22, 41, 42,
  '516TRFW021','516TRFW021','516TRFW021','616TRFW242','616TRFW242', &
!        50,51, 56, 62, 64,
  '516TRFW081','516TRFW081','516TRFW081','516TRFW060','516TRFW060', &
!        65,72, 73,81, 92,
  '516TRFW060','616TRFW073','616TRFW073','516TRFW081','616TRFW094', &
!        93,98,101,102,103,
  '616TRFW094','616TRFW094','516TRFW117','516TRFW108','516TRFW108', &
!       104,108,109,113,116,
  '516TRFW108','516TRFW108','516TRFW122','516TRFW108','516TRFW122', &
!       117,119,120,122,124,
  '516TRFW117','516TRFW117','516TRFW108','516TRFW122','516TRFW108', &
!       127,130,133,137,142,
  '516TRFW122','616TRFW264','516TRFW060','516TRFW122','616TRFW264', &
!       201,202,211,212,231,
  '516TRFW202','516TRFW202','616TRFW211','616TRFW211','616TRFW242', &
!       242,251,263,264,299,
  '616TRFW242','616TRFW242','616TRFW263','616TRFW264','616TRFW264', &
!       301,312,313,330,333,
  '500DVEW815','500DVEW312','500DVEW818','500DVEW801','500DVEW815', &
!       341,351,352,361,374,
  '500DVEW351','500DVEW351','500DVEW351','500DVEW361','500DVEW351', &
!       431,475,478,492,510,
  '500DVEW431','400DVEW475','500DVEW815','500DVEW351','616TRFW998', &
!       542,590,600,631,660,
  '500DVEW818','500DVEW351','500DVEW818','500DVEW631','500DVEW805', &
!       671,730,746,747,748,
  '500DVEW815','500DVEW805','616TRFW746','616TRFW747','616TRFW747', &
!       760,768,800,801,805,
  '500DVEW351','500DVEW815','500DVEW801','500DVEW801','500DVEW805', &
!       807,811,815,818,821,
  '500DVEW807','500DVEW811','500DVEW815','500DVEW818','500DVEW821', &
!       839,901,920,981,998,
  '500DVEW839','500DVEW815','500DVEW361','500DVEW981','500DVEW815', &
!       999
  '500DVEW815'/

!     GET EQUATION FROM EQNUME ARRAY
      DONE=0
      LAST = 91
      FIRST = 1
      DO 55, WHILE (DONE.EQ.0)
          HALF = (LAST - FIRST +1)/2 + FIRST
          IF(FIA(HALF) .EQ. SPEC)THEN
              DONE = HALF
          ELSEIF(FIRST .EQ. LAST) THEN
              ERRFLAG = 1
              DONE = -1
          ELSE IF (FIA(HALF) .LT. SPEC) THEN
              FIRST = HALF
          ELSE
              LAST = HALF - 1
          ENDIF
55 CONTINUE

      IF(DONE .LT. 0) THEN
          VOLEQ = EQNUM(LAST)
      ELSE
          VOLEQ = EQNUM(DONE)
      ENDIF

  RETURN
  END SUBROUTINE R5_PNWEQN

! FIA DEFAULT VOLUME EQUATION FOR WA AND OR
  SUBROUTINE R6_PNWEQN(FORST,SPEC,VAR,VOLEQ,ERRFLAG)
  CHARACTER*2 FORST,DIST,VAR
  CHARACTER*10 VOLEQ, EQNUMW(66), EQNUME(66)
  INTEGER ERRFLAG,SPEC,FIA(66)
  INTEGER DONE,HALF,FIRST,LAST
  DATA (FIA(I),I=1,66)/ &
      11, 15, 17, 19, 20, 21, 22, 41, 42, 64, &
      66, 72, 73, 81, 92, 93, 98,101,103,106, &
     108,113,116,117,119,122,130,202,211,212, &
     231,242,263,264,299,312,321,351,352,361, &
     375,376,431,475,478,492,500,510,542,590, &
     600,631,660,730,740,746,747,760,768,805, &
     815,818,920,981,998,999/
  DATA (EQNUMW(I),I=1,66)/ &
!         11, 15, 17, 19, 20,
  '632TRFW011','532TRFW015','632TRFW011','632TRFW011','532TRFW021', &
!         21, 22, 41, 42, 64,
  '532TRFW021','632TRFW011','532TRFW081','632TRFW242','532TRFW060', &
!         66, 72, 73, 81, 92,
  '532TRFW060','632TRFW073','632TRFW073','532TRFW081','632TRFW098', &
!         93, 98,101,103,106,
  '632TRFW098','632TRFW098','632TRFW108','632TRFW108','632TRFW264', &
!        108,113,116,117,119,
  '632TRFW108','632TRFW108','532TRFW122','532TRFW117','632TRFW108', &
!        122,130,202,211,212,
  '532TRFW122','632TRFW264','632TRFW202','632TRFW211','632TRFW211', &
!        231,242,263,264,299,
  '632TRFW242','632TRFW242','632TRFW263','632TRFW264','632TRFW264', &
!        312,321,351,352,361,
  '500DVEW312','500DVEW351','500DVEW351','500DVEW351','500DVEW361', &
!        375,376,431,475,478,
  '500DVEW351','500DVEW351','500DVEW431','400DVEW475','400DVEW475', &
!        492,500,510,542,590,
  '500DVEW351','500DVEW351','500DVEW351','500DVEW818','500DVEW351', &
!        600,631,660,730,740,
  '500DVEW351','500DVEW631','500DVEW351','500DVEW351','500DVEW351', &
!        746,747,760,768,805,
  '500DVEW351','500DVEW351','500DVEW351','500DVEW351','500DVEW805', &
!        815,818,920,981,998,
  '500DVEW815','500DVEW818','500DVEW351','500DVEW981','500DVEW351', &
!       999
  '500DVEW351'/
  DATA (EQNUME(I),I=1,66)/ &
!         11, 15, 17, 19, 20,21,
  '616TRFW019','616TRFW019','616TRFW019','616TRFW019','516TRFW021', &
!         21, 22, 41, 42, 64,
  '516TRFW021','616TRFW019','516TRFW081','616TRFW242','516TRFW060', &
!         66, 72, 73, 81, 92,
  '516TRFW060','616TRFW073','616TRFW073','516TRFW081','616TRFW098', &
!         93, 98,101,103,106,
  '616TRFW094','616TRFW094','616TRFW108','616TRFW108','616TRFW264', &
!        108,113,116,117,119,
  '616TRFW108','616TRFW108','616TRFW122','516TRFW117','616TRFW108', &
!        122,130,202,211,212,
  '616TRFW122','616TRFW264','616TRFW202','616TRFW211','616TRFW211', &
!        231,242,263,264,299,
  '616TRFW242','616TRFW242','616TRFW263','616TRFW264','616TRFW264', &
!        312,321,351,352,361,
  '500DVEW312','500DVEW351','500DVEW351','500DVEW351','500DVEW361', &
!        375,376,431,475,478,
  '500DVEW351','500DVEW351','500DVEW351','400DVEW475','400DVEW475', &
!        492,500,510,542,590,
  '500DVEW351','500DVEW351','500DVEW351','500DVEW818','500DVEW351', &
!        600,631,660,730,740,
  '500DVEW351','500DVEW631','500DVEW351','500DVEW351','500DVEW351', &
!        746,747,760,768,805,
  '500DVEW351','500DVEW351','500DVEW351','500DVEW351','500DVEW805', &
!        815,818,920,981,998,
  '500DVEW815','500DVEW818','500DVEW351','500DVEW981','500DVEW351', &
!       999
  '500DVEW351'/

  DONE=0
  LAST = 66
  FIRST = 1
!     Westside Variants
  IF(VAR.EQ.'PN' .OR. VAR.EQ.'WC' .OR. VAR.EQ.'NC' .OR. &
     VAR.EQ.'CA')THEN
!     GET EQUATION FROM EQNUMW ARRAY
      DO 65, WHILE (DONE.EQ.0)
          HALF = (LAST - FIRST +1)/2 + FIRST
          IF(FIA(HALF) .EQ. SPEC)THEN
              DONE = HALF
          ELSEIF(FIRST .EQ. LAST) THEN
              ERRFLAG = 1
              DONE = -1
          ELSE IF (FIA(HALF) .LT. SPEC) THEN
              FIRST = HALF
          ELSE
              LAST = HALF - 1
          ENDIF
65 CONTINUE

      IF(DONE .LT. 0) THEN
          VOLEQ = EQNUMW(LAST)
      ELSE
          VOLEQ = EQNUMW(DONE)
      ENDIF

!     FORESTS IN WA
!     NEED TO FIND FOREST NUMBER FOR WA!!!
    IF(FORST.EQ.'03' .OR. FORST.EQ.'05' .OR. FORST.EQ.'09' &
      .OR. FORST.EQ.'13' .OR. FORST.EQ.'17' .OR. FORST.EQ.'21')THEN
      IF(SPEC.EQ.72 .OR. SPEC.EQ.73)THEN
        VOLEQ='632TRFW073'
      ELSEIF(SPEC.EQ.122)THEN
        VOLEQ='632TRFW122'
      ELSEIF(SPEC.EQ.312 .OR. SPEC.EQ.361 .OR. SPEC.EQ.376 &
             .OR. SPEC.EQ.431 .OR. SPEC.EQ.492 &
             .OR. SPEC.EQ.542 .OR. SPEC.EQ.815)THEN
        VOLEQ='500DVEW351'
      ELSEIF(SPEC.EQ.351)THEN
        VOLEQ='616TRFW351'
      ENDIF
    ENDIF
!     EAST VARIANTS
  ELSE
!     GET EQUATION FROM EQNUME ARRAY
      DO 75, WHILE (DONE.EQ.0)
          HALF = (LAST - FIRST +1)/2 + FIRST
          IF(FIA(HALF) .EQ. SPEC)THEN
              DONE = HALF
          ELSEIF(FIRST .EQ. LAST) THEN
              ERRFLAG = 1
              DONE = -1
          ELSE IF (FIA(HALF) .LT. SPEC) THEN
              FIRST = HALF
          ELSE
              LAST = HALF - 1
          ENDIF
75 CONTINUE

      IF(DONE .LT. 0) THEN
          VOLEQ = EQNUME(LAST)
      ELSE
          VOLEQ = EQNUME(DONE)
      ENDIF

!     FORESTS IN WA
    IF(FORST.EQ.'03' .OR. FORST.EQ.'05' .OR. FORST.EQ.'09' &
      .OR. FORST.EQ.'13' .OR. FORST.EQ.'17' .OR. FORST.EQ.'21')THEN
      IF(SPEC.EQ.351)THEN
        VOLEQ='616TRFW351'
      ELSEIF(SPEC.EQ.312 .OR. SPEC.EQ.361 .OR. SPEC.EQ.542 &
              .OR. SPEC.EQ.815 .OR. SPEC.EQ.818)THEN
        VOLEQ='500DVEW351'
      ENDIF
    ENDIF
  ENDIF

  RETURN
  END SUBROUTINE R6_PNWEQN

