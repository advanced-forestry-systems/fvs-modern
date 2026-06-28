! test_gompsurv2.f90  -- unit test of the CORRECTED GOMPSURV (issue #75 fix).
! gsurv2 mirrors the fixed engine formula exactly: annual survival
! = 1 - exp(-exp(eta)), compounded over FINTL. Uses the affine cch map
! (cch_gompit = 0.062 + 0.0036*cch_hat) so inputs are on the fitted scale.
program test_gompsurv2
  implicit none
  integer, parameter :: MX = 600
  integer :: spcd(MX), nn, ios, ns, i
  real :: b0(MX), b1(MX), b2(MX), b3(MX), b4(MX), sben, sstr
  character(len=512) :: hdr
  open(10, file='greg_mortality_coefficients.csv', status='old')
  read(10,'(A)') hdr
  ns = 0
  do
    read(10,*,iostat=ios) spcd(ns+1), nn, b0(ns+1), b1(ns+1), b2(ns+1), b3(ns+1), b4(ns+1)
    if (ios /= 0) exit
    ns = ns + 1
    if (ns >= MX) exit
  end do
  close(10)
  open(20, file='surv2_fortran.csv')
  write(20,'(A)') 'SPCD,s_benign_ann,s_stress_ann'
  do i = 1, ns
    call gsurv2(b0(i),b1(i),b2(i),b3(i),b4(i), 0.6, cchg(10.0), 1.0, sben)
    call gsurv2(b0(i),b1(i),b2(i),b3(i),b4(i), 0.2, cchg(80.0), 1.0, sstr)
    write(20,'(I0,",",F12.9,",",F12.9)') spcd(i), sben, sstr
  end do
  close(20)
  write(*,'(A,I0)') 'wrote surv2_fortran.csv species=', ns
contains
  real function cchg(cch_hat)
    real cch_hat
    cchg = 0.062 + 0.0036*cch_hat
  end function cchg
end program test_gompsurv2

subroutine gsurv2(B0,B1,B2,B3,B4, CR, CCHV, FINTL, SURV)
  implicit none
  real B0,B1,B2,B3,B4, CR, CCHV, FINTL, SURV, CRC, CCHC, ETA, HZ, CTERM
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
  HZ = 1.0 - EXP(-EXP(ETA))                 ! annual gompit survival
  SURV = MAX(0.0, MIN(1.0, HZ)) ** FINTL    ! compound over the cycle
end subroutine gsurv2
