!== last modified  09-14-2007
  FUNCTION SF_CORR(JSP,geosub,TOTALH,HI,HJ)
!                     INTERNAL ROUTINES:
!                                          COR_C2
!                                          COR_OT
!                                          COR_BH
!                                          COR_WS
!
!     given 2 heights (hi and hj), estimate the correlation of the
!       errors in the corresponding dib's
  character*2 geosub
  INTEGER JSP
  REAL TOTALH,HI,HJ,SF_CORR,COR_C2,COR_BH,COR_OT,COR_WS,cor_ak

  if ( JSP.ge. 11  .and.  jsp.le. 21 ) then
      SF_CORR = COR_C2(JSP,geosub,TOTALH,HI, HJ)
      RETURN

  else if(JSP.eq.22) then
      SF_CORR=COR_BH(TOTALH,HI,HJ)
      RETURN
  else if(JSP.GE.23 .AND. JSP.LE.30) then
      SF_CORR=COR_OT(JSP,TOTALH,HI,HJ)
      RETURN
  else if(JSP.GE.31 .AND. JSP.LE.36) then
      SF_CORR=COR_AK(JSP,TOTALH,HI,HJ)
      RETURN
  else if( JSP.eq.1) then
!          SF_CORR=CORRH(HI,HJ)
      RETURN
  else  if ( JSP.eq.2) then
!          SF_CORR=CORRD(HI,HJ)
      RETURN
  else if (JSP.GE.3 .AND. JSP.LE.5) then
      SF_CORR=COR_WS(JSP,TOTALH,HI,HJ)
      RETURN
  endif

!      WRITE(*,11)
!11    format(' FUNCTION CORR GIVEN INVALID SPECIES.  STOP')
  RETURN
  END

