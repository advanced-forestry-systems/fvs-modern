! test_greghtdbh.f90 -- standalone Fortran evaluator for the CONUS native HT-DBH
! substitution (greghtdbh.f90 / GREGHTDBH). Transcribes the six Marshall model
! forms verbatim and, for every covered species and a DBH grid (1..50 in),
! writes predicted total height (feet) to htdbh_fortran.csv for cross-check
! against an R re-evaluation of each species' exact model.formula.
! NOTE: compiled real(kind=8) (double precision) to clear the 1e-6 validation
! gate; float32 only reaches ~4e-6 (cf. PR #95 fix/gregdghg-real8-precision).
program test_greghtdbh
  implicit none
  integer, parameter :: MX = 600
  integer :: spcd(MX), mdl(MX), nsp, ios, i, id
  real(kind=8) :: b1(MX), b2(MX), b3(MX)
  real(kind=8) :: dbh, ht
  character(len=2048) :: line
  integer :: imdl, ispcd
  real(kind=8) :: rb1, rb2, rb3
  open(10,file='conus_htdbh_coefficients.csv',status='old')
  read(10,'(A)') line
  nsp = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    call parse5(line, ispcd, imdl, rb1, rb2, rb3, ios)
    if (ios /= 0) cycle
    if (imdl < 1 .or. imdl > 6) cycle
    nsp = nsp + 1
    spcd(nsp)=ispcd; mdl(nsp)=imdl; b1(nsp)=rb1; b2(nsp)=rb2; b3(nsp)=rb3
    if (nsp >= MX) exit
  end do
  close(10)
  open(20,file='htdbh_fortran.csv')
  write(20,'(A)') 'SPCD,model_id,DBH,HT'
  do i=1,nsp
    do id=1,50
      dbh = dble(id)
      ht  = ghd_eval(mdl(i), b1(i), b2(i), b3(i), dbh)
      write(20,'(I0,",",I0,",",F6.1,",",ES22.14)') spcd(i), mdl(i), dbh, ht
    end do
  end do
  close(20)
  write(*,'(A,I0,A)') 'HT-DBH species=', nsp, ' -> htdbh_fortran.csv'
contains
  real(kind=8) function ghd_eval(m, p1, p2, p3, d)
    implicit none
    integer :: m
    real(kind=8) :: p1, p2, p3, d, bh
    bh = 4.5d0
    select case (m)
      case (1); ghd_eval = bh + exp(p1 + p2*d**(-1.0d0))
      case (2); ghd_eval = bh + exp(p1 + p2*d**p3)
      case (3); ghd_eval = bh + p1*(1.0d0-exp(p2*d))**p3
      case (4); ghd_eval = bh + exp(p1 + p2/(d+1.0d0))
      case (5); ghd_eval = bh + exp(p1 + p2/(d+p3))
      case (6); ghd_eval = bh + p1*(1.0d0-exp(p2*d**p3))
      case default; ghd_eval = -1.0d0
    end select
  end function ghd_eval
  subroutine parse5(ln, ifia, imod, v1, v2, v3, ierr)
    implicit none
    character(len=*) :: ln
    integer :: ifia, imod, ierr
    real(kind=8) :: v1, v2, v3
    integer :: p0, p1i, nf, lnn, ival
    real(kind=8) :: rval
    character(len=256) :: fld
    ifia=0; imod=0; v1=0d0; v2=0d0; v3=0d0; ierr=0
    lnn = len_trim(ln); p0 = 1; nf = 0
    do
      p1i = index(ln(p0:lnn), ',')
      if (p1i == 0) then
        fld = ln(p0:lnn)
      else
        fld = ln(p0:p0+p1i-2)
      end if
      nf = nf + 1
      if (len_trim(fld) == 0) then
        rval = 0d0; ival = 0
      else
        if (nf <= 2) then
          read(fld,*,iostat=ierr) ival
          if (ierr /= 0) return
        else
          read(fld,*,iostat=ierr) rval
          if (ierr /= 0) return
        end if
      end if
      select case (nf)
        case (1); ifia = ival
        case (2); imod = ival
        case (3); v1 = rval
        case (4); v2 = rval
        case (5); v3 = rval
      end select
      if (nf >= 5) exit
      if (p1i == 0) exit
      p0 = p0 + p1i
      if (p0 > lnn) exit
    end do
    if (nf < 4) ierr = 1
  end subroutine parse5
end program test_greghtdbh
