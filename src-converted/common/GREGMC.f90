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
!----------
INTEGER NGREGDG, NGREGHG
REAL    GDG(MAXSP,7), GHG(MAXSP,9), GMCW(MAXSP,3), GEMT, GTD, GELEV
LOGICAL LGREGDG, LGREGHG, GHAVE_DG(MAXSP), GHAVE_HG(MAXSP), GHAVE_MCW(MAXSP)
COMMON /GREGMR/ GDG, GHG, GMCW, GEMT, GTD, GELEV
COMMON /GREGMI/ NGREGDG, NGREGHG
COMMON /GREGML/ LGREGDG, LGREGHG, GHAVE_DG, GHAVE_HG, GHAVE_MCW
