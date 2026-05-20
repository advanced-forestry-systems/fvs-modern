C     Stub routines for uncompilable FVS volume functions.
C     These are safe no-ops: growth dynamics (DG, mortality, crown, H-D)
C     do not depend on volume. NSBE biomass is computed post-hoc.

      SUBROUTINE VOLS(ISESS,IFLAG)
C     Main volume computation entry point
      IMPLICIT NONE
      INTEGER ISESS,IFLAG
      RETURN
      END

      SUBROUTINE VOLEQDEF(VOLEQ,REGN,FORST,DIST,SPEC,PROD)
C     Volume equation definition lookup
      IMPLICIT NONE
      CHARACTER*(*) VOLEQ,PROD
      INTEGER REGN,FORST,DIST,SPEC
      VOLEQ = '          '
      RETURN
      END

      SUBROUTINE R9VOL(REGN,FORST,VOLEQ,MTEFLAG,NOLOGP,
     &  CUTFLG,BESSION,DBHOB,DRCOB,HTTOT,HTLOG,HT1PRD,
     &  HT2PRD,UPSHT1,UPSHT2,UPSD1,UPSD2,HTREF,AVGZ,
     &  HTUP,DPTS,VOL,LOGVOL,LOGDIA,LOGLEN,BOLHT,
     &  TLOGS,NUFLAG,PROD,LIVE,BA,SI,CTTYPE,HTTFLL,
     &  IDIST,BESSION2)
C     Region 9 volume computation
      IMPLICIT NONE
      INTEGER REGN,FORST,MTEFLAG,NOLOGP,CUTFLG,BESSION
      INTEGER NUFLAG,LIVE,IDIST,BESSION2
      CHARACTER*(*) VOLEQ,PROD,CTTYPE
      REAL DBHOB,DRCOB,HTTOT,HTLOG,HT1PRD,HT2PRD
      REAL UPSHT1,UPSHT2,UPSD1,UPSD2,HTREF,AVGZ,HTUP
      REAL BA,SI,HTTFLL,TLOGS
      REAL VOL(15),LOGVOL(7,20),LOGDIA(21,3)
      REAL LOGLEN(20),BOLHT(21),DPTS(*)
      INTEGER I,J
C     Initialize all volume outputs to zero
      DO I=1,15
        VOL(I) = 0.0
      ENDDO
      DO J=1,20
        LOGLEN(J) = 0.0
        DO I=1,7
          LOGVOL(I,J) = 0.0
        ENDDO
      ENDDO
      DO J=1,3
        DO I=1,21
          LOGDIA(I,J) = 0.0
        ENDDO
      ENDDO
      DO I=1,21
        BOLHT(I) = 0.0
      ENDDO
      TLOGS = 0.0
      NUFLAG = 0
      RETURN
      END

      SUBROUTINE VOLINITNVB(REGN,FORST,VOLEQ,DBHOB,HTTOT,
     &  VOL)
C     NVB volume initialization
      IMPLICIT NONE
      INTEGER REGN,FORST
      CHARACTER*(*) VOLEQ
      REAL DBHOB,HTTOT,VOL(15)
      INTEGER I
      DO I=1,15
        VOL(I) = 0.0
      ENDDO
      RETURN
      END

      SUBROUTINE VARVER(VERS)
C     Stub for variant version reporter
      IMPLICIT NONE
      CHARACTER*(*) VERS
      VERS = " "
      RETURN
      END
