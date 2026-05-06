!== last modified  09-14-2007
   SUBROUTINE SF_DFZ(JSP,GEOSUB,DBH,THT, H, DIBmod, Z, DIBact)
!      H       In    REAL*4    Section ht (0<= h < HT,  h not= BH)
!      DIBmod  In    Real*4    Base model prediction of dib (expectation)
!      Z       In    Real*4    Standardized offset
!      DIBact  Out   Real*4    Actual dib

  CHARACTER*2 GEOSUB
  INTEGER JSP
  REAL DBH,THT,H,DIBmod,Z,DIBact,SE_LNX
  REAL PROD_ZSE

  if(jsp.ge.11) then

     if (JSP .ge. 11 .and. jsp .le. 21) then
!                          Assume errors about dib are LOGNORMAL
        CALL VAR_C2(JSP,geosub,DBH,THT,H,SE_LNX)
     else if (JSP .ge.23 .and. jsp .le.30) then
!                          Assume errors about dib are LOGNORMAL
        CALL VAR_OT(JSP,DBH,THT,H,SE_LNX)
     else if(JSP.eq.22) then
!                              BLACK HILLS
        CALL VAR_BH(DBH,THT,H,SE_LNX)
     else if (JSP .ge.31 .and. jsp.le.36) then
!                          Assume errors about dib are LOGNORMAL
        CALL VAR_AK(JSP,DBH,THT,H,SE_LNX)
     endif
     PROD_ZSE = Z*SE_LNX
     if( abs(prod_ZSE) .gt. 5.0) then
        if (prod_ZSE .lt.  0.0) prod_ZSE=-5.0
        if (prod_ZSE .gt.  0.0) prod_ZSE= 5.0
     ENDIF
     DIBact = DIBmod * exp(prod_ZSE)
  ELSE
     if(JSP.eq. 1  .or. jsp.eq.2) then
!                          Assume errors about dib are NORMAL
!            CALL  VAR_STD( h, SE )
!            DIBact = DIBmod + Z * SE
     elseif(JSP.eq.3) then
!                          Special distributions (Johnson Sb)
        CALL SF_DFZ3(DBH,THT,H,DIBmod,DIBact,Z)
     else if(JSP.eq.4) then
        CALL SF_DFZ4(DBH,THT,H,DIBmod,DIBact,Z)
     else if(JSP.eq.5) then
        CALL SF_DFZ5(DBH,THT,H,DIBmod,DIBact,Z)
     else
!                           This is for unprogrammed species,
!                           lacking the three-point method.
        DIBact = DIBmod
      endif
  ENDIF

  RETURN
  END

