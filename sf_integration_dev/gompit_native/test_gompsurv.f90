! test_gompsurv.f90
! Faithful unit test of the engine GOMPSURV formula (base/gompmort.f90) against
! the generated coefficient CSV. The GSURV subroutine below is copied verbatim
! from GOMPSURV (same clamps, same eta, same hazard/survival), only swapping the
! GOMPMC common GB(ISPC,1:5) for passed coefficients. Reads
! greg_mortality_coefficients.csv exactly as GOMPLOAD does (skip header, then
! list-directed SPCD n b0 b1 b2 b3 b4) and writes surv for a grid of cr/cch.
program test_gompsurv
  implicit none
  integer, parameter :: MX = 600
  integer :: spcd(MX), nn, ios, ns, i, k
  real :: b0(MX), b1(MX), b2(MX), b3(MX), b4(MX)
  real :: crv(3), cchv(3), fint, surv
  character(len=512) :: header
  crv  = (/ 0.20, 0.50, 0.80 /)
  cchv = (/ 10.0, 50.0, 120.0 /)
  fint = 10.0
  open(10, file='greg_mortality_coefficients.csv', status='old')
  read(10,'(A)') header
  ns = 0
  do
    read(10,*,iostat=ios) spcd(ns+1), nn, b0(ns+1), b1(ns+1), b2(ns+1), b3(ns+1), b4(ns+1)
    if (ios /= 0) exit
    ns = ns + 1
    if (ns >= MX) exit
  end do
  close(10)
  open(20, file='surv_fortran.csv')
  write(20,'(A)') 'SPCD,cr,cch,fint,surv'
  do i = 1, ns
    do k = 1, 3
      call gsurv(b0(i),b1(i),b2(i),b3(i),b4(i), crv(k), cchv(k), fint, surv)
      write(20,'(I0,",",F6.3,",",F7.2,",",F5.1,",",F14.10)') spcd(i), crv(k), cchv(k), fint, surv
    end do
  end do
  close(20)
  write(*,'(A,I0,A)') 'wrote surv_fortran.csv for ', ns, ' species x 3 cr/cch points'
end program test_gompsurv

subroutine gsurv(B0,B1,B2,B3,B4, CR, CCHV, FINTL, SURV)
  implicit none
  real B0,B1,B2,B3,B4, CR, CCHV, FINTL, SURV
  real CRC, CCHC, ETA, HZ, CTERM
  CRC = CR
  if (CRC .lt. 1.0E-4) CRC = 1.0E-4
  if (CRC .gt. 1.0)    CRC = 1.0
  CCHC = CCHV
  if (CCHC .lt. 0.0) CCHC = 0.0
  if (CCHC .gt. 0.0) then
    CTERM = CCHC**B4
  else
    CTERM = 0.0
  end if
  ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM
  if (ETA .gt. 30.0)  ETA = 30.0
  if (ETA .lt. -30.0) ETA = -30.0
  HZ = EXP(ETA)
  SURV = EXP(-HZ*FINTL)
end subroutine gsurv
