# GREGDG / GREGHG native hook: build + validation record (2026-07-17)

Worktree: /fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-hook
Branch:   feat/greg-native-hook-20260717

## Wiring (all committed on this branch prior to this record)
- src-converted/common/GREGMC.f90    -- GREGMC common (LGREGDG/LGREGHG, GDG/GHG tables,
                                        GHAVE_DG/GHAVE_HG, per-stand GEMT/GTD/GELEV).
- src-converted/common/GREGKW.f90    -- keyword driver codes (IDGDRV/IHGDRV/ICRWDRV).
- src-converted/base/gregdghg.f90    -- GREGLOADDG (reads greg_dg_coefficients.csv,
                                        FIAJSP -> FVS species) + GREGDGV evaluator.
- src-converted/base/greghghg.f90    -- GREGLOADHG (HGDRIVER arms 1-4 or FVS_GREGHG env +
                                        FVS_GREGHG_COEF) + GREGHGV + crown-recession loaders.
- src-converted/ne/dgf.f90           -- DG hook: for LGREGDG.AND.GHAVE_DG(ISPC), replace the
                                        native WK2 with a 10-year GREGDGV loop, BRATIO to
                                        inside-bark DDS, WK2=ALOG(DDS) (COR omitted, per spec).
- src-converted/ne/htgf.f90          -- HG hook: for LGREGHG.AND.GHAVE_HG(ISPC), replace native
                                        HTG with a 10-year GREGHGV loop (CCFL from MCW-of-larger,
                                        CCH from GOMPCCH CCHT), bypassing BALMOD/SCALE/HTCON.
- src-converted/bin/FVSne_sourceList.txt -- includes base/gregdghg.f90 + base/greghghg.f90.

## Activation
- DG: env FVS_GREGDG=1 + FVS_GREGDG_COEF=config/greg_dg_coefficients.csv (9-col SPCD,n,B0..B6,
      elev+EMT form). Per-stand FVS_GREG_EMT / FVS_GREG_TD / FVS_GREG_ELEV scalars.
- HG: env FVS_GREGHG=1 + FVS_GREGHG_COEF=config/greg_hg_coefficients.csv (11-col SPCD,n,MX,B1..B8),
      or HGDRIVER keyword arm 1-4. Shares FVS_GREG_TD/EMT/ELEV.

## Build (sbatch job 12408842)
- build_greghook_ne.slurm -> bash deployment/scripts/build_fvs_executables.sh . ./bin-greghook ne
- Result: FVSne 7.9M, 557 objects, 98 skipped, exit 0.
- Symbol gate PASS: gregloaddg_, gregdgv_, gregloadhg_, greghgv_ all present (T).
- Binary: bin-greghook/FVSne (TEST path; production bin untouched).

## Validation
1. Unit harness (sf_integration_dev/gompit_native, test_gregdghg.f90 + validate_gregdghg.R),
   fixed state dbh=8,cr=0.5,ht=50,bal=80,ccfl=120,cch=0.4,elev=1500,td=25,emt=-15:
     DG  n=84  max|Fortran-R|=4.88e-12  median 0.136 in/yr  range [0.056, 0.375]
     HG  n=95  max|Fortran-R|=4.84e-11  median 1.322 ft/yr  range [0.106, 7.125]
     PASS match + PASS plausibility.
2. Native Fortran evaluator vs Python projector consumer (conus_greg_projector.py, branch
   feat/conus-greg-consumer-20260717), identical fixed state, annual increment:
     DG  n=84  max|native-projector|=2.37e-11
     HG  n=95  max|native-projector|=1.80e-09
   Confirms engine CSV coefficients == projector JSON block coefficients and identical forms.
3. NE integration smoke (net01.key, TEST binary):
   - Baseline (no env) vs Greg-activated both run to completion (10 cycles).
   - Activation confirmed in net01.out: "GREGDG enabled: 84 species" and "GREGHG enabled".
   - Trajectories diverge (2090 BA 85 -> 133, QMD 20.4 -> 34.4), i.e. the engine is consuming
     Greg's DG/HG rather than native NE-TWIGS.

## Remaining gaps
- Per-stand EMT/TD/ELEV are supplied as scalars via env (first-version path per spec). The
  full per-stand greg_emt_td_lookup (export to CSV keyed by stand id/fuzzed lat-lon) is not yet
  wired into GREGLOAD; stand-level A/B parity with the projector needs it.
- DG/HG cycle loops are hardcoded to 10 years (matches NE default FINT); not yet read from CONTRL.
- GOTCHA 1 (cch scale): HG feeds CCHT(I) from GOMPCCH; confirm it is the ORGANON 0-1 fraction
  (pre-affine), not the gompit affine CCHT, before production stand-level use.
- Stand-level native-vs-projector tolerance not closed (needs BAL/CCFL/CCH parity); the
  evaluator-level match (2.4e-11 / 1.8e-9) is the numeric ground truth here.

## Production gap closure (2026-07-17, later same day)

1. PER-STAND EMT/TD/ELEV LOOKUP (gap 1). New module src-converted/base/gregclim.f90
   (GREGCLIM_MOD: GREGCLIM_LOAD + GREGCLIM_APPLY). GREGLOADDG/GREGLOADHG now call
   GREGCLIM_LOAD once (reads FVS_GREG_CLIMATE_LOOKUP; CSV id,EMT,TD[,ELEV]; falls
   back to the FVS_GREG_EMT/TD/ELEV scalars, keeping that path intact). DGF and HTGF
   call GREGCLIM_APPLY per stand (match DBCN=STAND_CN then NPLT) to set GEMT/GTD/GELEV.
   Source: greg_emt_td_lookup.rds at fvs_stress/conus_eq_proj (1,861,183 stands,
   ClimateNA 1991-2020 EMT + TD=MWMT-MCMT). Exported per-stand CSV (STAND_CN,EMT,TD,ELEV,
   ELEV joined from standinit ELEVFT) to fvs_stress/conus_eq_proj/greg_climate_perstand.csv
   (n=1,861,183). Keyed by stand id only; no FIA coordinates released. Small NE sample:
   config/greg_climate_lookup_sample.csv; export script scripts/export_greg_climate.R.

2. CYCLE LENGTH (gap 2). DGF and HTGF Greg loops now iterate NCYR=NINT(FINT) annual
   sub-steps (FINT from PLOT common) instead of a hardcoded 10. HTGF drops the FINT/YR
   SCALE on the Greg branch (the true cycle length is now in NCYR). For a 10-year cycle
   NCYR=10, so 10-year results are byte-identical to the prior validated build.

3. CCHT CONVENTION (gap 3) -- CONFIRMED CORRECT, no code change. GOMPCCH fills
   CCHT(I)=CCH_A+CCH_B*cch_hat (0.062, 0.0036). The reference projector's cch_module.R
   stand_cch() returns pmax(CCH_A+CCH_B*raw,0) -- the SAME affine value -- and the
   projector feeds that single cch to BOTH hg_annual and gomp_surv_annual. The Fortran
   survival hook consumes CCHT (verified e.g. ec/morts.f90: CCHG=CCHT(I); CALL GOMPSURV).
   So the HG hook feeding CCHV=CCHT(I) is on the correct scale. The "pre-affine 0-1
   fraction" hypothesis was wrong: pre-affine cch_hat is the ORGANON crown-closure sum
   (0..~100, CCF-like), not 0-1; the affine map is exactly what puts it on the fit scale.

4. STAND-LEVEL PARITY (gap 4). TEST binary on net01.key (STAND S248112), Greg DG+HG
   active, per-stand climate from the sample CSV (EMT -12.5, TD 27.5, ELEV 2200 confirmed
   applied via GREGTRACE), FVS_GREG_TRACE=1. Per-tree engine GREGDGV/GREGHGV inputs
   (DBH,CR,HT,BAL / HT,CR,CCFL,CCH) fed to GregProjector.dg_annual/hg_annual:
     DG  n=52 unique  max|engine-projector| = 1.55e-7 in/yr (rel 3.6e-7)
     HG  n=26 unique  max|engine-projector| = 6.84e-7 ft/yr (rel 1.3e-7)
   Closes to float32 epsilon (engine GREGDGV is single precision); looser than the
   evaluator 2.4e-11/1.8e-9 exactly as expected from single-precision engine state.

Rebuild: sbatch 12409329, exit 0, 558 objects. Symbol gate PASS: gregloaddg_, gregdgv_,
gregloadhg_, greghgv_, plus __gregclim_mod_MOD_gregclim_load / _apply. Unit harness
re-run PASS: DG max|F-R|=4.88e-12, HG=4.84e-11.
