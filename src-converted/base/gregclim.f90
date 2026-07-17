!==============================================================================
!  gregclim.f90 -- GREGCLIM_MOD: per-stand EMT/TD/ELEV climate lookup for the
!  native Greg DG/HG hook (production gap 1).
!
!  GREGCLIM_LOAD reads a CSV keyed by stand id (STAND_CN / FVS stand id) with
!  EMT, TD and optional ELEV columns, once, from FVS_GREG_CLIMATE_LOOKUP. It
!  also captures the scalar FVS_GREG_EMT / FVS_GREG_TD / FVS_GREG_ELEV env
!  values as the fallback, so the pre-existing scalar activation path keeps
!  working unchanged when no lookup file is supplied.
!
!  GREGCLIM_APPLY resolves the current stand: it matches the stand id string
!  against DBCN (STAND_CN) first, then NPLT (FVS stand id), and on a hit
!  overrides EMT/TD (and ELEV when the CSV carried it). Called once per stand
!  from DGF and HTGF, ahead of the species loop, so GEMT/GTD/GELEV in /GREGMR/
!  and /GREGHR/ carry the correct per-stand climate into GREGDGV / GREGHGV.
!
!  No FIA true coordinates are used or stored here: the table is keyed by
!  stand id only. The upstream export (extract_emt_td.R) samples ClimateNA
!  normals at the fuzzed standinit lat/lon and emits STAND_CN,EMT,TD[,ELEV].
!
!  CSV format: a header line, then rows "ID,EMT,TD[,ELEV]". Comma OR whitespace
!  delimited; ELEV column optional (falls back to the scalar / native ELEV).
!==============================================================================
MODULE GREGCLIM_MOD
  IMPLICIT NONE
  INTEGER, SAVE :: NGCLIM = 0
  LOGICAL, SAVE :: LGCLIM = .FALSE.        ! usable per-stand lookup loaded
  LOGICAL, SAVE :: LGCLIM_TRIED = .FALSE.  ! load attempted (guard: once)
  LOGICAL, SAVE :: GCLIM_HASELEV = .FALSE. ! CSV carried an ELEV column
  REAL,    SAVE :: GCLIM_EMT0 = 0.0, GCLIM_TD0 = 0.0, GCLIM_ELEV0 = 0.0
  CHARACTER(LEN=40), ALLOCATABLE, SAVE :: GCLIM_ID(:)
  REAL, ALLOCATABLE, SAVE :: GCLIM_EMT(:), GCLIM_TD(:), GCLIM_ELEV(:)
CONTAINS

  SUBROUTINE GREGCLIM_LOAD(JOUT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: JOUT
    CHARACTER(LEN=256)  :: CPATH, CVAL
    CHARACTER(LEN=1024) :: LINE
    INTEGER :: IOS, U, N, K, IC
    REAL    :: E, T, EL
    CHARACTER(LEN=40) :: IDS
    IF (LGCLIM_TRIED) RETURN
    LGCLIM_TRIED = .TRUE.
    ! scalar env fallbacks (same vars the scalar path already reads)
    CALL GETENV('FVS_GREG_EMT',  CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GCLIM_EMT0
    CALL GETENV('FVS_GREG_TD',   CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GCLIM_TD0
    CALL GETENV('FVS_GREG_ELEV', CVAL); IF (CVAL.NE.' ') READ(CVAL,*,IOSTAT=IOS) GCLIM_ELEV0
    CALL GETENV('FVS_GREG_CLIMATE_LOOKUP', CPATH)
    IF (CPATH.EQ.' ') RETURN                 ! no lookup -> scalar path only
    U = 72
    OPEN(UNIT=U, FILE=CPATH, STATUS='OLD', IOSTAT=IOS)
    IF (IOS.NE.0) THEN
      WRITE(JOUT,*) 'GREGCLIM: cannot open ',TRIM(CPATH),'; scalar climate only.'
      RETURN
    ENDIF
    ! pass 1: count non-blank data rows (skip header)
    N = 0
    READ(U,'(A)',IOSTAT=IOS) LINE
    DO
      READ(U,'(A)',IOSTAT=IOS) LINE
      IF (IOS.NE.0) EXIT
      IF (LEN_TRIM(LINE).EQ.0) CYCLE
      N = N + 1
    ENDDO
    IF (N.LE.0) THEN
      CLOSE(U); WRITE(JOUT,*) 'GREGCLIM: empty lookup ',TRIM(CPATH); RETURN
    ENDIF
    ALLOCATE(GCLIM_ID(N), GCLIM_EMT(N), GCLIM_TD(N), GCLIM_ELEV(N))
    GCLIM_ELEV = 0.0
    ! pass 2: parse (comma OR whitespace; 3 or 4 fields, ELEV optional)
    REWIND(U)
    READ(U,'(A)',IOSTAT=IOS) LINE
    K = 0
    DO
      READ(U,'(A)',IOSTAT=IOS) LINE
      IF (IOS.NE.0) EXIT
      IF (LEN_TRIM(LINE).EQ.0) CYCLE
      DO IC=1,LEN(LINE)
        IF (LINE(IC:IC).EQ.',') LINE(IC:IC)=' '
      ENDDO
      EL = -9999.0
      READ(LINE,*,IOSTAT=IOS) IDS, E, T, EL
      IF (IOS.NE.0) THEN
        EL = -9999.0
        READ(LINE,*,IOSTAT=IOS) IDS, E, T
      ENDIF
      IF (IOS.NE.0) CYCLE
      K = K + 1
      GCLIM_ID(K)  = ADJUSTL(IDS)
      GCLIM_EMT(K) = E
      GCLIM_TD(K)  = T
      IF (EL.GT.-9998.0) THEN
        GCLIM_ELEV(K) = EL
        GCLIM_HASELEV = .TRUE.
      ENDIF
    ENDDO
    CLOSE(U)
    NGCLIM = K
    LGCLIM = (K.GT.0)
    WRITE(JOUT,*) 'GREGCLIM: per-stand climate lookup loaded, n=',K, &
                  ' hasELEV=',GCLIM_HASELEV,' from ',TRIM(CPATH)
    RETURN
  END SUBROUTINE GREGCLIM_LOAD

  SUBROUTINE GREGCLIM_APPLY(IDCN, IDPLT, EMT, TD, ELEV, LHIT)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: IDCN, IDPLT
    REAL, INTENT(INOUT) :: EMT, TD, ELEV
    LOGICAL, INTENT(OUT) :: LHIT
    INTEGER :: K
    CHARACTER(LEN=40) :: A, B
    LHIT = .FALSE.
    IF (.NOT.LGCLIM) RETURN
    A = ADJUSTL(IDCN)
    B = ADJUSTL(IDPLT)
    DO K=1,NGCLIM
      IF (GCLIM_ID(K).EQ.A .OR. GCLIM_ID(K).EQ.B) THEN
        EMT = GCLIM_EMT(K)
        TD  = GCLIM_TD(K)
        IF (GCLIM_HASELEV) ELEV = GCLIM_ELEV(K)
        LHIT = .TRUE.
        RETURN
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE GREGCLIM_APPLY

END MODULE GREGCLIM_MOD
