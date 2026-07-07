program diag
  implicit none
  integer, parameter :: MX=600
  integer :: spcd(MX), nn, ios, ns, k, id, ic, ih
  real :: b0(MX),b1(MX),b2(MX),b3(MX),b4(MX)
  real :: dbhg(6), crg(3), cchg(3), surv, s10
  character(len=512) :: header
  integer :: focus(3)
  dbhg=(/0.6,1.0,1.5,2.0,5.0,10.0/)
  crg =(/0.20,0.45,0.80/)
  cchg=(/0.0,50.0,120.0/)
  focus=(/12,94,241/)
  open(10,file="greg_mortality_coefficients.csv",status="old")
  read(10,"(A)") header
  ns=0
  do
    read(10,*,iostat=ios) spcd(ns+1),nn,b0(ns+1),b1(ns+1),b2(ns+1),b3(ns+1),b4(ns+1)
    if(ios/=0) exit
    ns=ns+1; if(ns>=MX) exit
  end do
  close(10)
  write(*,*) "SPCD dbh cr cch surv_ann surv_10yr mort10_pct"
  do id=1,ns
    do k=1,3
      if(spcd(id)==focus(k)) then
        do ih=1,6
         do ic=1,3
          call gsurv(b0(id),b1(id),b2(id),b3(id),b4(id),dbhg(ih),crg(ic),50.0,surv,s10)
          write(*,900) spcd(id),dbhg(ih),crg(ic),50.0,surv,s10,(1.0-s10)*100.0
         end do
        end do
      end if
    end do
  end do
900 format(I5,2X,F5.2,2X,F5.2,2X,F6.1,2X,F9.6,2X,F9.6,2X,F7.2)
end program diag

subroutine gsurv(B0,B1,B2,B3,B4,DBHV,CR,CCHV,SURV,S10)
  implicit none
  real B0,B1,B2,B3,B4,DBHV,CR,CCHV,SURV,S10
  real CRC,CCHC,ETA,HZ,CTERM,DBHC
  CRC=CR; if(CRC.lt.1.0E-4)CRC=1.0E-4; if(CRC.gt.1.0)CRC=1.0
  CCHC=CCHV; if(CCHC.lt.0.0)CCHC=0.0
  if(CCHC.gt.0.0)then; CTERM=CCHC**B4; else; CTERM=0.0; endif
  DBHC=DBHV; if(DBHC.lt.1.0)DBHC=1.0
  ETA=B0+B1*(CRC+0.01)**B2+B3*CTERM
  if(ETA.gt.30.0)ETA=30.0; if(ETA.lt.-30.0)ETA=-30.0
  HZ=1.0-EXP(-EXP(ETA))
  SURV=MAX(0.0,MIN(1.0,HZ))
  S10=SURV**10.0
end subroutine gsurv
