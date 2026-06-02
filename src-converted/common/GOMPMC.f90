!----------
!  GOMPMC -- Greg Johnson gompit mortality, in-engine state.
!
!  Shared state for the gompit mortality option (annual hazard on crown ratio
!  and crown closure at tree tip). Must be INCLUDEd after PRGPRM.f90 so MAXSP
!  and MAXTRE are defined.
!
!    LGOMP  -- master switch; .TRUE. when gompit mortality is active (set by
!              GOMPLOAD from the FVS_GOMPIT environment variable).
!    NGOMP  -- number of FVS species with a fitted gompit row.
!    GB     -- per-FVS-species coefficient table, GB(ISPC,1:5) = b0,b1,b2,b3,b4
!              resolved from the FIA-keyed coefficient file via FIAJSP.
!    GHAVE  -- .TRUE. for FVS species that carry a fitted gompit row; gompit
!              owns their mortality. Unfit species keep native background.
!    CCHT   -- per-tree crown closure at tip for the current cycle, on the
!              gompit cch scale (affine-mapped). Filled by GOMPCCH each cycle.
!----------
INTEGER NGOMP, GGRP(MAXSP)
REAL GB(MAXSP,5), CCHT(MAXTRE)
LOGICAL LGOMP, GHAVE(MAXSP)
COMMON /GOMPMR/ GB, CCHT
COMMON /GOMPMI/ NGOMP, GGRP
COMMON /GOMPML/ LGOMP, GHAVE
