!----------
!  GOMPMC -- Greg Johnson gompit mortality, in-engine state.
!
!  Shared state for the gompit mortality option (annual hazard on crown ratio
!  and crown closure at tree tip). Must be INCLUDEd after PRGPRM.f90 so MAXSP
!  and MAXTRE are defined.
!
!    LGOMP  -- master switch; .TRUE. when gompit mortality is active (set by
!              GOMPLOAD from the FVS_GOMPIT env var or the GOMPMORT keyword).
!    LGOMPKW-- set .TRUE. by the GOMPMORT keyword (via GOMPON in initre);
!              GOMPLOAD then activates in MORCON. Initialised .FALSE. in
!              BLOCK DATA GOMPBD.
!    NGOMP  -- number of FVS species with a fitted gompit row.
!    GB     -- per-FVS-species coefficient table, GB(ISPC,1:5) = b0,b1,b2,b3,b4
!              resolved from the FIA-keyed coefficient file via FIAJSP.
!    GHAVE  -- .TRUE. for FVS species that carry a fitted gompit row; gompit
!              owns their mortality. Unfit species keep native background.
!    CCHT   -- per-tree crown closure at tip for the current cycle, on the
!              gompit cch scale (affine-mapped). Filled by GOMPCCH each cycle.
!    GDBHMIN-- small-DBH guard threshold (inches). Trees with 0<=DBH<GDBHMIN
!              get the young-cohort survival floor instead of raw gompit, so
!              seedling cohorts do not collapse. Trees with DBH>=GDBHMIN are
!              BYTE-IDENTICAL to prior gompit. Default 1.0; override with the
!              FVS_GOMP_DBHMIN env var. Set <=0 to disable the guard entirely.
!    GSFLOOR-- annual-survival floor applied below GDBHMIN (fraction, 0..1).
!              Default 0.95 (~5%/yr, ~40%/decade max mortality). Override with
!              the FVS_GOMP_SFLOOR env var. The floor is a MAX on annual
!              survival, so gompit still governs where it predicts survival
!              above the floor.
!----------
INTEGER NGOMP, GGRP(MAXSP)
REAL GB(MAXSP,5), CCHT(MAXTRE), GDBHMIN, GSFLOOR
LOGICAL LGOMP, LGOMPKW, GHAVE(MAXSP)
COMMON /GOMPMR/ GB, CCHT, GDBHMIN, GSFLOOR
COMMON /GOMPMI/ NGOMP, GGRP
COMMON /GOMPML/ LGOMP, LGOMPKW, GHAVE
