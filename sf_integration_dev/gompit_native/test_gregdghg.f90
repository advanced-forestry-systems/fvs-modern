! test_gregdghg.f90 -- standalone Fortran evaluators for Greg's deployed DG and HG
! increments in Greg's PC climate form, transcribed from conus_eq_projector_greg_par.R
! (dg_annual, hg_annual). Reads the emitted coefficient CSVs and writes per-species increments
! at one representative tree state, to cross-check against the R projector forms.
! NOTE: recompiled real(kind=8) (double precision) to clear the 1e-6 validation
! gate; kind=8 chosen for consistency with existing double-precision declarations
! elsewhere in the engine (e.g. src-converted/base/apisubs.f90).
program test_gregdghg
  implicit none
  integer, parameter :: MX = 600
  integer :: dgsp(MX), hgsp(MX), ndg, nhg, nn, ios, i
  real(kind=8) :: d0(MX),d1(MX),d2(MX),d3(MX),d4(MX),d5(MX),d6(MX)
  real(kind=8) :: h0(MX),h1(MX),h2(MX),h3(MX),h4(MX),h5(MX),h6(MX),h7(MX)
  real(kind=8) :: dbh,cr,ht,bal,ccfl,cch,dd0,td,ppt,dd18,pc1,pc2, dg, hg
  character(len=800) :: hdr
  ! representative mid-size tree state
  dbh=8.0d0; cr=0.5d0; ht=50.0d0; bal=80.0d0; ccfl=120.0d0; cch=0.4d0
  dd0=1600.0d0; td=28.0d0; ppt=230.0d0; dd18=3500.0d0
  pc1=0.252d0*dd0+0.002d0*td-0.035d0*ppt+0.967d0*dd18
  pc2=0.882d0*dd0+0.015d0*td+0.420d0*ppt-0.215d0*dd18
  open(10,file='greg_dg_coefficients.csv',status='old'); read(10,'(A)') hdr; ndg=0
  do
    read(10,*,iostat=ios) dgsp(ndg+1),nn,d0(ndg+1),d1(ndg+1),d2(ndg+1),d3(ndg+1),d4(ndg+1),d5(ndg+1),d6(ndg+1)
    if(ios/=0) exit; ndg=ndg+1; if(ndg>=MX) exit
  end do
  close(10)
  open(11,file='greg_hg_coefficients.csv',status='old'); read(11,'(A)') hdr; nhg=0
  do
    read(11,*,iostat=ios) hgsp(nhg+1),nn,h0(nhg+1),h1(nhg+1),h2(nhg+1),h3(nhg+1),h4(nhg+1),h5(nhg+1),h6(nhg+1),h7(nhg+1)
    if(ios/=0) exit; nhg=nhg+1; if(nhg>=MX) exit
  end do
  close(11)
  open(20,file='dg_fortran.csv'); write(20,'(A)') 'SPCD,dg'
  do i=1,ndg
    call dgcalc(d0(i),d1(i),d2(i),d3(i),d4(i),d5(i),d6(i), dbh,cr,ht,bal,pc1,pc2, dg)
    write(20,'(I0,",",ES18.10)') dgsp(i), dg
  end do
  close(20)
  open(21,file='hg_fortran.csv'); write(21,'(A)') 'SPCD,hg'
  do i=1,nhg
    call hgcalc(h0(i),h1(i),h2(i),h3(i),h4(i),h5(i),h6(i),h7(i), ht,cr,ccfl,cch,pc1,pc2, hg)
    write(21,'(I0,",",ES18.10)') hgsp(i), hg
  end do
  close(21)
  write(*,'(A,I0,A,I0)') 'DG species=', ndg, '  HG species=', nhg
end program test_gregdghg

subroutine dgcalc(B0,B1,B2,B3,B4,B5,B6, dbh,cr,ht,bal,pc1,pc2, dg)
  implicit none
  real(kind=8) :: B0,B1,B2,B3,B4,B5,B6, dbh,cr,ht,bal,pc1,pc2, dg, z
  z = B0 + B1*log((dbh+1.0d0)**2/(cr*ht+1.0d0)**B3) + B2*bal**B4/log(dbh+2.7d0) + B5*pc1 + B6*pc2
  if (z .gt. 5.0d0)   z = 5.0d0
  if (z .lt. -30.0d0) z = -30.0d0
  dg = exp(z); if (dg .lt. 0.0d0) dg = 0.0d0
end subroutine dgcalc

subroutine hgcalc(mx,b1,b2,b3,b4,b5,b6,b7, ht,cr,ccfl,cch,pc1,pc2, hg)
  implicit none
  real(kind=8) :: mx,b1,b2,b3,b4,b5,b6,b7, ht,cr,ccfl,cch,pc1,pc2, hg, crp,cchp,arg
  crp = cr;  if (crp .lt. 1.0d-4) crp = 1.0d-4
  cchp = cch; if (cchp .lt. 0.0d0) cchp = 0.0d0
  arg = -b1*ht - b4*ccfl - b5*cchp**0.5d0 + b6*pc1 + b7*pc2
  hg = mx*b1*b2*crp**b3*exp(arg)*(1.0d0-exp(-b1*ht))**(b2-1.0d0)
  if (hg .lt. 0.0d0) hg = 0.0d0
end subroutine hgcalc
