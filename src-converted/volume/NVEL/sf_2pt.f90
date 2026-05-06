!== last modified  4-9-2002
  SUBROUTINE SF_2PT(JSP,GEOSUB,SETOPT,DBH,TOTALH,DBT_BH, &
                                 F,RHFW,RFLW,TAPCOE)
!                                                               SF_2PT
!                                  separate entry:              SF_2PTX1
!                                  separate entry:              SF_2PTX
!     Stem form   Specify the 2-point input: DBH and Total height
!     DBH_arg   input   R*4    DBH
!     HT_arg    input   R*4    Total tree height
!     ASP       input   A2     User-defined species code
!                                  (blank is OK  IF specieas already set)
!     IERR      output  I       0  No Error
!                               1  Unrecognized species code
!                               2  Invalid DBH
!                               3  Invalid height
!                               4  This species not available in metric.
!                               5  Supplemental arguments can not be used
!                                  because regression coef file not
!                                  available for this species.
!                               6  Region code not recognized.
  CHARACTER*2 GEOSUB
!      LOGICAL IVSTAT_2, IVSTAT_S
!      LOGICAL METRIC, SUPPLEMENT
!      REAL*4 SUPV_ARG(NSUP)
  INTEGER  SETOPT(6),IERR
!      LOGICAL  SF_MCOE(NSP)
  REAL*4 RHFW(4),RFLW(6),TAPCOE(12),F,YHAT_BH
  REAL DBH, TOTALH, DBH_IB, DBT_BH, BH, RH, slope
  REAL DBT_USER,SF_YHAT
  INTEGER JSP, INEEDSL

  DBT_USER = dbt_bh
  bh = 4.5
!      DO 10 I=1,NSUP
!         Supv(I) =supv_arg(I)
!         if (supv(i).ne. 0.)  SUPPLEMENT=.TRUE.
!10    Continue

  INEEDSL=0

  IF( JSP.eq.0)then
      IERR=1
      RETURN
  ENDIF

!      if (JSP.ge.3 .and. METRIC) then
!        IERR=4
!        RETURN
!      ENDIF

!      IVSTAT_2 = .false.
!      IVSTAT_S = .false.
  F=1.0
!       find shape parameters
  CALL SF_SHP(JSP,GEOSUB,SETOPT,DBH,TOTALH,RFLW,RHFW,DBT_BH, &
                DBT_USER,DBH_IB)

!          use formulas by Flewelling and Raynes to  calc. taper coef's

  CALL SF_TAPER (RHFW,RFLW,TAPCOE)

!          calculate scaling factor f
  RH = BH/TOTALH

  YHAT_BH = SF_YHAT(JSP,RH,totalh,ineedsl,slope,RHFW,RFLW,TAPCOE,F)
  F  = DBH_IB / YHAT_BH

!                  if segment volumes may be used, define extra sections
!           IF(DBH_IB .gt. 6.5) then
!           IVEXTRA=3
!           DVEXTRA(1) = 6.0
!           DVEXTRA(2) = 4.0
!           DVEXTRA(3) = 0.0
!           if (metric)  then
!                 dvextra(1)= dvextra(1)*2.54
!                 dvextra(2)= dvextra(2)*2.54
!               ENDIF
!           HVEXTRA(1)=0.
!           HVEXTRA(2)=0.
!           HVEXTRA(3)=RHI1*TOTALH
!         else
!           IVEXTRA=0.0
!         ENDIF

  RETURN
  END
