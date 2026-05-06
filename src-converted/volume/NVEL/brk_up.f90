!== last modified  06-24-2005
  SUBROUTINE BRK_UP(JSP,GEOSUB,DBHOB,THT,DBTBH,HTUP,DIB,DOB,DBT)
!             given a section height (HTUP) and dib, calculate the ratio:
!                            DBT/DOB
  character*2 geosub
  INTEGER JSP
  REAL DBHOB,THT,HTUP,DIB,DBT,DBTBH,BARK_R,RATIO,DOB
  REAL BRK_WS,BRK_UPA2

  IF (JSP.GE.22 .AND. JSP.LE.30) THEN
       DOB = DIB
      CALL BRK_OT(JSP,geosub,DBHOB,DOB,HTUP,DBTBH,DIB,dbt)
  elseif (jsp.ge.11 .and. jsp.lt.22) then
    IF(DIB.le. 0.) then
      dbt=0.0
      DOB = 0
    ELSE
      RATIO = BRK_UPA2(JSP,DBHOB,THT,DBTBH,HTUP,dib)
      if(ratio.le.0.0 .or. ratio .ge. 1.0) then
         dbt=0.0
         DOB = DIB
      else
         dbt = RATIO/(1.0-RATIO)*DIB
         DOB = DIB+DBT
      endif
    ENDIF

!      ELSEIF (JSP.GE.31 .AND. JSP.LE.34) THEN
!          CALL BRK_AK(JSP,DBHOB,DOB,HTUP,DBTBH,DIB,DBT)

  else if(JSP.GE.3 .AND. JSP.LE.5) THEN
!          write(*,*)' bark',DBHOB,dbtbh,HTUP
     BARK_R = BRK_WS(JSP,DBHOB,THT,DBTBH,HTUP)
     dbt = dib * (bark_r / (1.0-bark_r))
     dob = dib + dbt
  ENDIF

!      if(jsp.eq.1) then
!        brk_upa=bark_upH(HTUP)
!      else  if(jsp.eq.2) then
!        brk_upa=bark_upD(HTUP)
!      endif
!       write(6,11)
!11     format(' Subroutine Bark_up given invalid species.  STOP')

  RETURN
  end



