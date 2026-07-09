!----------
!  GREGMC -- Greg Johnson native diameter-growth (and future HG) substitution state.
!  INCLUDE after PRGPRM.f90 so MAXSP/MAXTRE are defined. Mirrors GOMPMC.f90.
!
!    LGREGDG  -- master switch; .TRUE. when Greg DG substitution is active (set by
!                GREGLOADDG from the FVS_GREGDG env var).
!    NGREGDG  -- number of FVS species with a fitted Greg DG row.
!    GDG      -- per-FVS-species DG coefficient table, GDG(ISPC,1:7) = B0..B6,
!                resolved from the FIA-keyed coefficient file via FIAJSP.
!    GHAVE_DG -- .TRUE. for FVS species that carry a fitted Greg DG row; Greg DG
!                owns their diameter increment. Unfit species keep native NE-TWIGS.
!    GEMT/GTD/GELEV -- per-stand climate (extreme min temp, temp diff) and elevation
!                (feet), resolved once per stand by GREGLOADDG.
!    GDGMAX   -- per-FVS-species maximum diameter (in) for the size-based DG
!                deceleration in GREGDGV. Read from an optional 10th column
!                (DBHMAX) of the coefficient file; if the column is absent the
!                array is left at a large sentinel so decel==1 (old behaviour).
!----------
INTEGER NGREGDG
REAL    GDG(MAXSP,7), GDGMAX(MAXSP), GDD0, GTD, GPPT_SM, GDD18, GPC1, GPC2
LOGICAL LGREGDG, GHAVE_DG(MAXSP)
COMMON /GREGMR/ GDG, GDGMAX, GDD0, GTD, GPPT_SM, GDD18, GPC1, GPC2
COMMON /GREGMI/ NGREGDG
COMMON /GREGML/ LGREGDG, GHAVE_DG
