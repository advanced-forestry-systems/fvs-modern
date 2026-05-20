! varver_stub.f90
!
! Stub for VARVER (variant version reporting routine). FVS variant code
! (e.g. acd/dgdriv.f90, acd/morts.f90) calls VARVER to retrieve a variant-
! version identifier, but no variant in upstream USDA FVS actually defines
! VARVER. The shared-library build tolerates this because .so allows
! unresolved symbols at link time, but standalone executable builds reject
! it. This stub provides the missing definition.
!
! Originally bundled inside econ/econ_stubs.f90 in the holoros fork. Split
! out 2026-05-16 so that the econ stubs can be conditionally excluded for
! variants with the real econ implementation (ACD, NE, CS, LS, SN) while
! VARVER remains universally available.
!
! When a future upstream commit lands an actual variant-version routine,
! this stub should be deleted and the build script entry removed.

subroutine VARVER(ESSION)
  integer ESSION
  return
end subroutine VARVER
