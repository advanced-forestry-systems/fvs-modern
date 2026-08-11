! Invariance + guard unit test. GSURV_OLD = pre-guard GOMPSURV verbatim.
! GSURV_NEW = post-guard GOMPSURV verbatim (DBH-gated floor). Prove:
!  (1) for DBH>=GDBHMIN: NEW == OLD bit-for-bit  (established invariance)
!  (2) for DBH<GDBHMIN: NEW = max(OLD, floor)     (seedling guard)
program invtest
  implicit none
  integer, parameter :: MX=600
  integer :: spcd(MX),nn,ios,ns,id,ic,ih
  real :: b0(MX),b1(MX),b2(MX),b3(MX),b4(MX)
  real :: dbhg(5),crg(3),cch,so,sn,fint,GDBHMIN,GSFLOOR
  character(len=512) :: header
  integer :: nviol_est, nfloor, ntot
  dbhg=(/0.6,1.0,1.5,5.0,15.0/); crg=(/0.20,0.45,0.80/); cch=0.30; fint=10.0
  GDBHMIN=1.0; GSFLOOR=0.95
  open(10,file="greg_mortality_coefficients.csv",status="old"); read(10,"(A)") header
  ns=0
  do; read(10,*,iostat=ios) spcd(ns+1),nn,b0(ns+1),b1(ns+1),b2(ns+1),b3(ns+1),b4(ns+1)
    if(ios/=0)exit; ns=ns+1; if(ns>=MX)exit; end do
  close(10)
  nviol_est=0; nfloor=0; ntot=0
  do id=1,ns
    do ih=1,5
      do ic=1,3
        call gold(b0(id),b1(id),b2(id),b3(id),b4(id),crg(ic),cch,fint,so)
        call gnew(b0(id),b1(id),b2(id),b3(id),b4(id),dbhg(ih),crg(ic),cch,fint,GDBHMIN,GSFLOOR,sn)
        ntot=ntot+1
        if(dbhg(ih).ge.GDBHMIN) then
          if(sn.ne.so) nviol_est=nviol_est+1   ! must be exactly equal
        else
          if(sn.gt.so+1e-9) nfloor=nfloor+1     ! floor lifted survival
          if(sn.lt.so-1e-6) nviol_est=nviol_est+1 ! floor must never lower survival
        endif
      end do
    end do
  end do
  write(*,*) "species=",ns," cases=",ntot
  write(*,*) "established-tree (DBH>=",GDBHMIN,") NEW/=OLD violations =",nviol_est," (want 0)"
  write(*,*) "seedling cases where floor lifted survival =",nfloor
  ! show one worst-case seedling species before/after
  do id=1,ns
    if(spcd(id)==241) then
      call gold(b0(id),b1(id),b2(id),b3(id),b4(id),0.20,cch,fint,so)
      call gnew(b0(id),b1(id),b2(id),b3(id),b4(id),0.6,0.20,cch,fint,GDBHMIN,GSFLOOR,sn)
      write(*,*) "SPCD241 seedling(DBH0.6,cr0.2): OLD surv10=",so," NEW surv10=",sn
      write(*,*) "  OLD mort10=",(1.0-so)*100.,"%  NEW mort10=",(1.0-sn)*100.,"%"
      call gnew(b0(id),b1(id),b2(id),b3(id),b4(id),15.0,0.20,cch,fint,GDBHMIN,GSFLOOR,sn)
      write(*,*) "SPCD241 established(DBH15,cr0.2): OLD=",so," NEW=",sn," (equal? ",so.eq.sn,")"
    endif
  end do
end program

subroutine gold(B0,B1,B2,B3,B4,CR,CCHV,FINTL,SURV)
  implicit none; real B0,B1,B2,B3,B4,CR,CCHV,FINTL,SURV,CRC,CCHC,ETA,HZ,CTERM
  CRC=CR; if(CRC.lt.1e-4)CRC=1e-4; if(CRC.gt.1.0)CRC=1.0
  CCHC=CCHV; if(CCHC.lt.0.0)CCHC=0.0
  if(CCHC.gt.0.0)then; CTERM=CCHC**B4; else; CTERM=0.0; endif
  ETA=B0+B1*(CRC+0.01)**B2+B3*CTERM
  if(ETA.gt.30.0)ETA=30.0; if(ETA.lt.-30.0)ETA=-30.0
  HZ=1.0-EXP(-EXP(ETA)); SURV=MAX(0.0,MIN(1.0,HZ))**FINTL
end subroutine
subroutine gnew(B0,B1,B2,B3,B4,DBHV,CR,CCHV,FINTL,GDBHMIN,GSFLOOR,SURV)
  implicit none
  real B0,B1,B2,B3,B4,DBHV,CR,CCHV,FINTL,GDBHMIN,GSFLOOR,SURV
  real CRC,CCHC,ETA,HZ,CTERM,SFLR
  CRC=CR; if(CRC.lt.1e-4)CRC=1e-4; if(CRC.gt.1.0)CRC=1.0
  CCHC=CCHV; if(CCHC.lt.0.0)CCHC=0.0
  if(CCHC.gt.0.0)then; CTERM=CCHC**B4; else; CTERM=0.0; endif
  ETA=B0+B1*(CRC+0.01)**B2+B3*CTERM
  if(ETA.gt.30.0)ETA=30.0; if(ETA.lt.-30.0)ETA=-30.0
  HZ=1.0-EXP(-EXP(ETA)); SURV=MAX(0.0,MIN(1.0,HZ))**FINTL
  if(GDBHMIN.gt.0.0 .and. DBHV.ge.0.0 .and. DBHV.lt.GDBHMIN)then
    SFLR=GSFLOOR**FINTL; if(SURV.lt.SFLR)SURV=SFLR
  endif
end subroutine
