! test_gregdghg.f90 -- standalone Fortran evaluators for Greg's deployed DG and HG
! increments, transcribed verbatim from conus_eq_projector_greg.R (dg_annual,
! hg_annual). Reads the emitted coefficient CSVs and writes per-species increments
! at one representative tree state, to cross-check against the R projector forms.
program test_gregdghg
  implicit none
  integer, parameter :: MX = 600
  integer :: dgsp(MX), hgsp(MX), ndg, nhg, nn, ios, i
  real :: d0(MX),d1(MX),d2(MX),d3(MX),d4(MX),d5(MX),d6(MX)
  real :: h0(MX),h1(MX),h2(MX),h3(MX),h4(MX),h5(MX),h6(MX),h7(MX),h8(MX)
  real :: dbh,cr,ht,bal,ccfl,cch,elev,td,emt, dg, hg
  character(len=800) :: hdr
  ! representative mid-size tree state
  dbh=8.0; cr=0.5; ht=50.0; bal=80.0; ccfl=120.0; cch=0.4; elev=1500.0; td=25.0; emt=-15.0
  open(10,file='greg_dg_coefficients.csv',status='old'); read(10,'(A)') hdr; ndg=0
  do
    read(10,*,iostat=ios) dgsp(ndg+1),nn,d0(ndg+1),d1(ndg+1),d2(ndg+1),d3(ndg+1),d4(ndg+1),d5(ndg+1),d6(ndg+1)
    if(ios/=0) exit; ndg=ndg+1; if(ndg>=MX) exit
  end do
  close(10)
  open(11,file='greg_hg_coefficients.csv',status='old'); read(11,'(A)') hdr; nhg=0
  do
    read(11,*,iostat=ios) hgsp(nhg+1),nn,h0(nhg+1),h1(nhg+1),h2(nhg+1),h3(nhg+1),h4(nhg+1),h5(nhg+1),h6(nhg+1),h7(nhg+1),h8(nhg+1)
    if(ios/=0) exit; nhg=nhg+1; if(nhg>=MX) exit
  end do
  close(11)
  open(20,file='dg_fortran.csv'); write(20,'(A)') 'SPCD,dg'
  do i=1,ndg
    call dgcalc(d0(i),d1(i),d2(i),d3(i),d4(i),d5(i),d6(i), dbh,cr,ht,bal,elev,emt, dg)
    write(20,'(I0,",",ES18.10)') dgsp(i), dg
  end do
  close(20)
  open(21,file='hg_fortran.csv'); write(21,'(A)') 'SPCD,hg'
  do i=1,nhg
    call hgcalc(h0(i),h1(i),h2(i),h3(i),h4(i),h5(i),h6(i),h7(i),h8(i), ht,cr,ccfl,cch,elev,td,emt, hg)
    write(21,'(I0,",",ES18.10)') hgsp(i), hg
  end do
  close(21)
  write(*,'(A,I0,A,I0)') 'DG species=', ndg, '  HG species=', nhg
end program test_gregdghg

subroutine dgcalc(B0,B1,B2,B3,B4,B5,B6, dbh,cr,ht,bal,elev,emt, dg)
  implicit none
  real B0,B1,B2,B3,B4,B5,B6, dbh,cr,ht,bal,elev,emt, dg, z
  z = B0 + B1*log((dbh+1.0)**2/(cr*ht+1.0)**B3) + B2*bal**B4/log(dbh+2.7) + B5*elev + B6*emt
  if (z .gt. 5.0)   z = 5.0
  if (z .lt. -30.0) z = -30.0
  dg = exp(z); if (dg .lt. 0.0) dg = 0.0
end subroutine dgcalc

subroutine hgcalc(mx,b1,b2,b3,b4,b5,b6,b7,b8, ht,cr,ccfl,cch,elev,td,emt, hg)
  implicit none
  real mx,b1,b2,b3,b4,b5,b6,b7,b8, ht,cr,ccfl,cch,elev,td,emt, hg, crp,cchp,tdp,arg
  crp = cr;  if (crp .lt. 1.0e-4) crp = 1.0e-4
  cchp = cch; if (cchp .lt. 0.0) cchp = 0.0
  tdp = td;  if (tdp .lt. 0.0) tdp = 0.0
  arg = -b1*ht - b4*ccfl - b8*cchp**0.5 - b5*elev + b6*sqrt(tdp) + b7*emt
  hg = mx*b1*b2*crp**b3*exp(arg)*(1.0-exp(-b1*ht))**(b2-1.0)
  if (hg .lt. 0.0) hg = 0.0
end subroutine hgcalc
