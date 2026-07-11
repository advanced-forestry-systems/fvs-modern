!----------
!  GREGMC -- Greg Johnson DG, HG, and crown-recession substitution state.
!  INCLUDE after PRGPRM.f90 so MAXSP is defined. Mirrors GOMPMC.f90.
!
!  ==== DG variables (PC-form, 9-column CSV B0..B6 + DBHMAX) ====
!
!    LGREGDG  -- master switch; .TRUE. when Greg DG substitution is active (set by
!                GREGLOADDG from the FVS_GREGDG env var or DGDRIVER keyword).
!    NGREGDG  -- number of FVS species with a fitted Greg DG row.
!    GDG      -- per-FVS-species DG coefficient table, GDG(ISPC,1:7) = B0..B6,
!                resolved from the FIA-keyed coefficient file via FIAJSP.
!    GDGMAX   -- per-species maximum diameter (in) for size-based DG deceleration.
!    GHAVE_DG -- .TRUE. for FVS species that carry a fitted Greg DG row.
!    GDD0/GTD/GPPT_SM/GDD18 -- per-stand ClimateNA inputs for PC-form DG.
!    GPC1/GPC2 -- principal component climate scalars derived from the above.
!
!  ==== HG variables (Chapman-Richards ORGANON form, 11-column CSV) ====
!
!    LGREGHG   -- master switch for Greg HG substitution.
!    NGREGHG   -- number of FVS species with a fitted Greg HG row.
!    GHG       -- GHG(ISPC,1:9): (1)=MX asymptote, (2:9)=B1..B8.
!    GHAVE_HG  -- .TRUE. for species that carry a fitted Greg HG row.
!    GMCW      -- GMCW(ISPC,1:4): (1)=FORM code, (2:4)=A,B,C for MCW quadratic.
!    GHAVE_MCW -- .TRUE. for species that carry a fitted MCW row.
!    GEMT/GTD/GELEV -- per-stand climate for HG (EMT deg C, TD deg C, elevation ft).
!                      GTD is shared with DG (same FVS_GREG_TD env var).
!
!  ==== Crown recession variables (3-param per-species model) ====
!
!    LGREG_CRW  -- master switch for Greg crown recession substitution.
!    NGREG_CRW  -- number of FVS species with a fitted crown recession row.
!    GCRW       -- GCRW(ISPC,1:3) = b0, b1, b2.
!    GHAVE_CRW  -- .TRUE. for species that carry a fitted crown recession row.
!----------
! -- DG --
INTEGER NGREGDG
REAL    GDG(MAXSP,7), GDGMAX(MAXSP), GDD0, GTD, GPPT_SM, GDD18, GPC1, GPC2
LOGICAL LGREGDG, GHAVE_DG(MAXSP)
COMMON /GREGMR/ GDG, GDGMAX, GDD0, GTD, GPPT_SM, GDD18, GPC1, GPC2
COMMON /GREGMI/ NGREGDG
COMMON /GREGML/ LGREGDG, GHAVE_DG
! -- HG --
INTEGER NGREGHG
REAL    GHG(MAXSP,9), GMCW(MAXSP,4), GEMT, GELEV
LOGICAL LGREGHG, GHAVE_HG(MAXSP), GHAVE_MCW(MAXSP)
COMMON /GREGHR/ GHG, GMCW, GEMT, GELEV
COMMON /GREGHC/ NGREGHG
COMMON /GREGHL/ LGREGHG, GHAVE_HG, GHAVE_MCW
! -- Crown recession --
INTEGER NGREG_CRW
REAL    GCRW(MAXSP,3)
LOGICAL LGREG_CRW, GHAVE_CRW(MAXSP)
COMMON /GREGCR/ GCRW
COMMON /GREGCC/ NGREG_CRW
COMMON /GREGCL/ LGREG_CRW, GHAVE_CRW
