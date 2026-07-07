!----------
!  GREGHDMC -- Marshall CONUS native height-diameter (HT-DBH) substitution state.
!  INCLUDE after PRGPRM.f90 so MAXSP is defined. Mirrors GREGMC.f90 (Greg DG state).
!
!    LGREGHTDBH -- master switch; .TRUE. when the CONUS HT-DBH substitution is
!                  active (set by GREGLOADHD from the FVS_GREGHTDBH env var).
!    NGREGHD    -- number of FVS species carrying a fitted CONUS HT-DBH row.
!    GHDMOD     -- per-FVS-species model_id (1..6) selecting one of the six
!                  Marshall HT-DBH forms; 0 for uncovered species.
!    GHD        -- per-FVS-species HT-DBH coefficient table, GHD(ISPC,1:3)=B1,B2,B3.
!                  B3 is unused (0) for models 1 and 4.
!    GHAVE_HD   -- .TRUE. for FVS species that carry a fitted CONUS HT-DBH row; the
!                  native Wykoff/Curtis-Arney prediction is replaced for these.
!----------
INTEGER NGREGHD, GHDMOD(MAXSP)
REAL    GHD(MAXSP,3)
LOGICAL LGREGHTDBH, GHAVE_HD(MAXSP)
COMMON /GREGHDR/ GHD
COMMON /GREGHDI/ NGREGHD, GHDMOD
COMMON /GREGHDL/ LGREGHTDBH, GHAVE_HD
