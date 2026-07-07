! Find seedling collapse cases (OLD surv10 < 0.10) across a cch sweep that a
! young DENSE cohort can produce, and show guard rescue. Focus NE/SN conifers.
program collapse
  implicit none
  integer, parameter :: MX=600
  integer :: spcd(MX),nn,ios,ns,id,ic,ih
  real :: b0(MX),b1(MX),b2(MX),b3(MX),b4(MX)
  real :: crg(3),cchsw(6),so,sn,fint,GDBHMIN,GSFLOOR
  character(len=512) :: header
  crg=(/0.10,0.30,0.60/); cchsw=(/0.062,0.5,2.0,5.0,15.0,40.0/); fint=10.0
  GDBHMIN=1.0; GSFLOOR=0.95
  open(10,file="greg_mortality_coefficients.csv",status="old"); read(10,"(A)") header
  ns=0
  do; read(10,*,iostat=ios) spcd(ns+1),nn,b0(ns+1),b1(ns+1),b2(ns+1),b3(ns+1),b4(ns+1)
    if(ios/=0)exit; ns=ns+1; if(ns>=MX)exit; end do
  close(10)
  write(*,*) "SPCD cr cch  OLDsurv10 OLDmort% -> NEWsurv10 NEWmort%  (DBH=0.6 seedling)"
  do id=1,ns
    if(spcd(id)==12 .or. spcd(id)==94 .or. spcd(id)==95 .or. spcd(id)==97 .or. spcd(id)==241) then
      do ic=1,3
        do ih=1,6
          call gold(b0(id),b1(id),b2(id),b3(id),b4(id),crg(ic),cchsw(ih),fint,so)
          call gnew(b0(id),b1(id),b2(id),b3(id),b4(id),0.6,crg(ic),cchsw(ih),fint,GDBHMIN,GSFLOOR,sn)
          if(so.lt.0.5) write(*,900) spcd(id),crg(ic),cchsw(ih),so,(1.-so)*100.,sn,(1.-sn)*100.
        end do
      end do
    endif
  end do
900 format(I5,1X,F4.2,1X,F6.2,1X,F9.5,1X,F7.2,1X,"->",1X,F9.5,1X,F7.2)
end program
subroutine gold(B0,B1,B2,B3,B4,CR,CCHV,FINTL,SURV)
  implicit none; real B0,B1,B2,B3,B4,CR,CCHV,FINTL,SURV,CRC,CCHC,ETA,HZ,CTERM
  CRC=CR; if(CRC.lt.1e-4)CRC=1e-4; if(CRC.gt.1.0)CRC=1.0
  CCHC=CCHV; if(CCHC.lt.0.0)CCHC=0.0
  if(CCHC.gt.0.0)then; CTERM=CCHC**B4; else; CTERM=0.0; endif
  ETA=B0+B1*(CRC+0.01)**B2+B3*CTERM; if(ETA.gt.30.)ETA=30.; if(ETA.lt.-30.)ETA=-30.
  HZ=1.0-EXP(-EXP(ETA)); SURV=MAX(0.0,MIN(1.0,HZ))**FINTL
end subroutine
subroutine gnew(B0,B1,B2,B3,B4,DBHV,CR,CCHV,FINTL,GDBHMIN,GSFLOOR,SURV)
  implicit none; real B0,B1,B2,B3,B4,DBHV,CR,CCHV,FINTL,GDBHMIN,GSFLOOR,SURV
  real CRC,CCHC,ETA,HZ,CTERM,SFLR
  CRC=CR; if(CRC.lt.1e-4)CRC=1e-4; if(CRC.gt.1.0)CRC=1.0
  CCHC=CCHV; if(CCHC.lt.0.0)CCHC=0.0
  if(CCHC.gt.0.0)then; CTERM=CCHC**B4; else; CTERM=0.0; endif
  ETA=B0+B1*(CRC+0.01)**B2+B3*CTERM; if(ETA.gt.30.)ETA=30.; if(ETA.lt.-30.)ETA=-30.
  HZ=1.0-EXP(-EXP(ETA)); SURV=MAX(0.0,MIN(1.0,HZ))**FINTL
  if(GDBHMIN.gt.0.0 .and. DBHV.ge.0.0 .and. DBHV.lt.GDBHMIN)then
    SFLR=GSFLOOR**FINTL; if(SURV.lt.SFLR)SURV=SFLR; endif
end subroutine
