SUBROUTINE DBSCASE(IFORSURE)
IMPLICIT NONE
!
! DBSQLITE $Id$
!
!
!     PURPOSE: TO POPULATE A DATABASE WITH THE PROGNOSIS MODEL
!              OUTPUT.
!
!     INPUT: IFORSURE - 1 NEED CONNECTION, 0 CHECK IF CONNECTION IS
!                         NEEDED.
!                       2 UPDATE SAMPLINGWT AND GROUPS IN CASES TABLE
!
!---
!OMMONS
!
!
INCLUDE 'PRGPRM.f90'
!
!
INCLUDE 'CONTRL.f90'
!
!
INCLUDE 'DBSCOM.f90'
!
!
INCLUDE 'KEYCOM.f90'
!
!
INCLUDE 'OPCOM.f90'
!
!
INCLUDE 'OUTCOM.f90'
!
!
INCLUDE 'PLOT.f90'
!
!
INCLUDE 'INCLUDESVN.f90'
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_ADDCOLIFABSENT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_BIND_DOUBLE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_BIND_INT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_BIND_TEXT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_CLOSE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLCNT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLDOUBLE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLINT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLISNULL &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLNAME &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLREAL &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLTEXT &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_COLTYPE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_ERRMSG &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_EXEC &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_FINALIZE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_OPEN &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_PREPARE &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_RESET &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_STEP &
   !$ ATTRIBUTES DLLIMPORT :: FSQL3_TABLEEXISTS &
   !_WIN64) &
   !$ ATTRIBUTES ALIAS:'_FSQL3_ADDCOLIFABSENT' :: FSQL3_ADDCOLIFABSENT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_BIND_DOUBLE'    :: FSQL3_BIND_DOUBLE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_BIND_INT'       :: FSQL3_BIND_INT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_BIND_TEXT'      :: FSQL3_BIND_TEXT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_CLOSE'          :: FSQL3_CLOSE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLCNT'         :: FSQL3_COLCNT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLDOUBLE'      :: FSQL3_COLDOUBLE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLINT'         :: FSQL3_COLINT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLISNULL'      :: FSQL3_COLISNULL &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLNAME'        :: FSQL3_COLNAME &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLREAL'        :: FSQL3_COLREAL &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLTEXT'        :: FSQL3_COLTEXT &
   !$ ATTRIBUTES ALIAS:'_FSQL3_COLTYPE'        :: FSQL3_COLTYPE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_ERRMSG'         :: FSQL3_ERRMSG &
   !$ ATTRIBUTES ALIAS:'_FSQL3_EXEC'           :: FSQL3_EXEC &
   !$ ATTRIBUTES ALIAS:'_FSQL3_FINALIZE'       :: FSQL3_FINALIZE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_OPEN'           :: FSQL3_OPEN &
   !$ ATTRIBUTES ALIAS:'_FSQL3_PREPARE'        :: FSQL3_PREPARE &
   !$ ATTRIBUTES ALIAS:'_FSQL3_RESET'          :: FSQL3_RESET &
   !$ ATTRIBUTES ALIAS:'_FSQL3_STEP'           :: FSQL3_STEP &
   !$ ATTRIBUTES ALIAS:'_FSQL3_TABLEEXISTS'    :: FSQL3_TABLEEXISTS

!OMMONS
!---
INTEGER fsql3_tableexists,fsql3_exec,fsql3_errmsg
INTEGER,parameter:: MaxStringLen=255

CHARACTER*2000 SQLStmtStr
CHARACTER*10  DATO, REV
CHARACTER*8   TIM, SVN
CHARACTER*(*) CFN
CHARACTER(len=MaxStringLen) TIMESTAMP
INTEGER IFORSURE, IFORSR, I, KODE, IRCODE
CHARACTER(len=36) CID
!-----
!     If output DB is active, indicated by presence of a nonblank CASEID,
!     update the SamplingWt and Groups in cases table that has been
!     determined since creation of initial record.
!     If IFORSURE is 2 and CASEID is still blank at this time, DB output
!     is not active (no case record to update) and just return.
!     The call is from INITRE in the PROCESS keyword section.
!-----
!      write(*,*) 'in DBSCASE TOP - IFORESURE, CASEID',IFORSURE,CASEID

IF(IFORSURE.EQ.2 .AND. CASEID.NE.' ') THEN
  IF (LENSLS.EQ.-1) THEN
    SLSET =""
  ENDIF
  WRITE(SQLStmtStr,*) "UPDATE FVS_Cases SET SamplingWt = ",SAMWT, &
      ", Groups = '", TRIM(ADJUSTL(SLSET)),"'", &
      ", StandID = '", TRIM(ADJUSTL(NPLT)),"'", &
      " WHERE CaseID = '",CASEID,"';"

!        write(*,*) 'in DBSCASE - ',SQLStmtStr

  IRCODE = fsql3_exec(IoutDBref,trim(SQLStmtStr)//CHAR(0))

  IF (IRCODE.ne.0) then
    IRCODE = fsql3_errmsg(IoutDBref,SQLStmtStr,LEN(SQLStmtStr)-1)
    PRINT *," IoutDBref=",IinDBref, &
               " ErrMsg =",SQLStmtStr(:IRCODE)
  ENDIF

  RETURN
ELSEIF(IFORSURE.EQ.2 .AND. CASEID.EQ.' ') THEN
  RETURN
ENDIF

!-----
!     CHECK TO SEE IF WE ARE NEEDING TO CONTINUE
!-----
IFORSR=IFORSURE
IF(IFORSR.EQ.0) THEN
  IF(ISUMARY.GE.1.OR. &
        ICOMPUTE.GE.1.OR. &
        ICALIB.GE.1.OR. &
        IATRTLIST.GE.1.OR. &
        ITREELIST.GE.1.OR. &
        ICUTLIST.GE.1.OR. &
        IDM1.GE.1.OR.IDM2.GE.1.OR.IDM3.GE.1.OR. &
        IDM5.GE.1.OR.IDM6.GE.1.OR. &
        IPOTFIRE.GE.1.OR. &
        IPOTFIREC.GE.1.OR. &
        IFUELS.GE.1.OR. &
        ICMRPT.GE.1.OR. &
        ICHRPT.GE.1.OR. &
        ISTRCLAS.GE.1.OR. &
        IFUELC.GE.1.OR. &
        IBURN.GE.1.OR. &
        IMORTF.GE.1.OR. &
        ISSUM.GE.1.OR. &
        ISDET.GE.1.OR. &
        ICANPR.GE.1.OR. &
        IDWDVOL.GE.1.OR. &
        IDWDCOV.GE.1.OR. &
        IBMMAIN.GE.1.OR. &
        IBMBKP.GE.1.OR. &
        IBMTREE.GE.1.OR. &
        IBMVOL.GE.1.OR. &
        ISTATS1.GE.1.OR.ISTATS2.GE.1.OR. &
        IREG1.GE.1.OR.IREG2.GE.1.OR.IREG3.GE.1.OR. &
        IREG4.GE.1.OR.IREG5.GE.1) IFORSR = 1
 ENDIF
 IF(IFORSR.EQ.0) RETURN
!---------
!     IF ALREADY HAVE A CURRENT CASE NUMBER, JUST BAIL
!---------
IF (CASEID.NE."") RETURN
KODE=0
!---------
!     MAKE SURE WE HAVE AN OPEN CONNECTION
!---------
IF (IoutDBref.EQ.-1) CALL DBSOPEN(.TRUE.,.FALSE.,KODE)
IF (KODE.EQ.1) THEN
  ICOMPUTE  = 0
  ICALIB    = 0
  ISUMARY   = 0
  IATRTLIST = 0
  ITREELIST = 0
  ICUTLIST  = 0
  IDM1      = 0
  IDM2      = 0
  IDM3      = 0
  IDM5      = 0
  IDM6      = 0
  IPOTFIRE  = 0
  IPOTFIREC = 0
  IFUELS    = 0
  ICMRPT    = 0
  ICHRPT    = 0
  ISTRCLAS  = 0
  IFUELC    = 0
  IBURN     = 0
  IMORTF    = 0
  ISSUM     = 0
  ISDET     = 0
  ICANPR    = 0
  IDWDVOL   = 0
  IDWDCOV   = 0
  IBMMAIN   = 0
  IBMBKP    = 0
  IBMTREE   = 0
  IBMVOL    = 0
  ISTATS1   = 0
  ISTATS2   = 0
  RETURN
ENDIF

IRCODE = fsql3_tableexists(IoutDBref,"FVS_Cases"//CHAR(0))
IF(IRCODE.EQ.0) THEN
  SQLStmtStr="CREATE TABLE FVS_Cases"// &
                 " (CaseID text primary key,"// &
                 "Stand_CN text not null,"// &
                 "StandID text not null,"// &
                 "MgmtID text,"// &
                 "RunTitle text,"// &
                 "KeywordFile text,"// &
                 "SamplingWt real,"// &
                 "Variant text,"// &
                 "Version text,"// &
                 "RV text,"// &
                 "Groups text,"// &
                 "RunDateTime text)"//CHAR(0)
  IRCODE = fsql3_exec(IoutDBref,SQLStmtStr)
  IF (IRCODE .NE. 0) RETURN
ENDIF
!---------
!     CREATE ENTRY FROM DATA FOR FVSRUN TABLE
!---------
CALL UUIDGEN(CASEID)
!     CREATE DATETIME

CALL GRDTIM(DATO,TIM)
TIMESTAMP = DATO(7:10)//'-'//DATO(1:5)//'-'//TIM

!     GET THE REVISION.
CALL REVISE(VARACD,REV)

!     STRIP 3-CHARACTER EXTENSION (IF PRESENT) FROM KEYFNAME

KEYFNAME = TRIM(KEYFNAME)
I = INDEX(".KEY",KEYFNAME)
IF (I.GT.0 .AND. I.EQ.LEN_TRIM(KEYFNAME)-3) THEN
  KEYFNAME = KEYFNAME(1:I-4)
ENDIF
IF (LENSLS.EQ.-1) THEN
  SLSET =""
ELSE
  CALL REMOVEQUOTES(SLSET)
ENDIF
CALL REMOVEQUOTES(DBCN)
CALL REMOVEQUOTES(NPLT)
CALL REMOVEQUOTES(MGMID)
CALL REMOVEQUOTES(ITITLE)
CALL REMOVEQUOTES(KEYFNAME)
IF (KEYFNAME.EQ.' ') THEN
  KEYFNAME='Unknown'
ELSE
  CALL REMOVEQUOTES(KEYFNAME)
ENDIF

WRITE(SQLStmtStr,*)"INSERT INTO FVS_Cases", &
    " (CaseID,Stand_CN,StandID,MgmtID,RunTitle,KeywordFile,", &
    "SamplingWt,Variant,Version,RV,Groups,RunDateTime) ", &
    "VALUES('",CASEID,"','", &
    TRIM(ADJUSTL(DBCN)),"','", &
    TRIM(ADJUSTL(NPLT)),"','", &
    TRIM(ADJUSTL(MGMID)),"','", &
    TRIM(ADJUSTL(ITITLE)),"','", &
    TRIM(ADJUSTL(KEYFNAME)),"',", &
    SAMWT,",'",VARACD,"','", &
    TRIM(ADJUSTL(SVN)),"','", &
    TRIM(ADJUSTL(REV)),"','", &
    TRIM(ADJUSTL(SLSET)),"','", &
    TRIM(ADJUSTL(TIMESTAMP)),"');"

!     write(*,*) 'DBSCASE_SQL: ',trim(SQLStmtStr)

IRCODE = fsql3_exec(IoutDBref,trim(SQLStmtStr)//CHAR(0))
IF (IRCODE.ne.0) then
    IRCODE = fsql3_errmsg(IoutDBref,SQLStmtStr,LEN(SQLStmtStr)-1)
    PRINT *," IoutDBref=",IinDBref, &
                " ErrMsg =",SQLStmtStr(:IRCODE)
ENDIF
RETURN
!
!     CALLED BY FILOPN: ENTRY TO SAVE THE KEYWORD FILE NAME AND TO SET
!                       THE DEFAULT DBS CONNECTIONS
!
ENTRY DBSVKFN (CFN)
I=LEN_TRIM(CFN)
IF (I.GT.LEN(KEYFNAME)) THEN
   KEYFNAME=CFN(1:4)//'...'//CFN(I+8-LEN(KEYFNAME):)
ELSE
   KEYFNAME = CFN
ENDIF
RETURN
!
!======================================================================
!     ENTRY FOR WWPBM, FETCHING CASEID
!    (NOTE: THE WWPBM NEEDS TO KNOW AND SAVE (INTERNALLY) CASEID, BECAUSE
!     IT IS DOING ITS DB-WRITING FROM WITHIN ITS OWN INTERNAL STAND LOOP
!
ENTRY DBSWW2(CID)
CID=CASEID
RETURN
!======================================================================
END

SUBROUTINE REMOVEQUOTES (CS)
IMPLICIT NONE
CHARACTER(LEN=*) CS
INTEGER I
DO I=1,LEN(CS)
  IF (CS(I:I) .EQ. "'" .OR. ICHAR(CS(I:I)) .LT. 32) CS(I:I)=' '
ENDDO
RETURN
END
