# GREGDG / GREGHG: native engine hooks for Greg's diameter and height growth (option 1)

Design for wiring Greg's deployed DG and HG into the FVS growth loop, mirroring the working GOMPMORT pattern
(base/gompmort.f90 + GOMPMC common + env/keyword activation + a one-line substitution in the host
subroutine). The evaluators below are already unit-validated against the projector (DG to 8e-8, HG to 4e-6;
see sf_integration_dev/gompit_native/test_gregdghg.f90 + validate_gregdghg.R).

## Forms (validated; from conus_eq_projector_greg.R)
DG annual diameter increment (in/yr), coeffs B0..B6:
    z  = B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7) + B5*elev + B6*EMT
    dg = exp(clamp(z,-30,5)), >= 0
HG annual height increment (ft/yr), coeffs B0=max_ht, b1..b8:
    dht = mx*b1*b2*cr^b3*exp(-b1*ht -b4*ccfl -b8*cch^0.5 -b5*elev +b6*sqrt(TD) +b7*EMT)*(1-exp(-b1*ht))^(b2-1)

## Engine substitution points (NE; generalizes to all eastern variants)
- DGF(DIAM) computes DDS (change in squared diameter) per tree into WK2. Greg gives an annual dbh increment g,
  so for a covered tree: dnew = dbh + g*FINT; DDS = dnew*dnew - dbh*dbh; store into WK2 in place of the native
  NE-TWIGS DDS. FINT = cycle length in years (CONTRL).
- HTGF computes the periodic height increment into HTG. Greg gives an annual height increment dht, so for a
  covered tree: HTG(i) = dht*FINT, replacing the native HTCALC increment.
Both substitutions are guarded per species by a GHAVE flag, exactly like morts.f90 does
(IF(LGREG.AND.GHAVE_DG(ISPC)) ... ), so unfit species keep the native equations.

## GREGMC common block (new; pattern of GOMPMC.f90)
    INTEGER NGDG, NGHG
    REAL    GDG(MAXSP,7)        ! B0..B6 per FVS species
    REAL    GHG(MAXSP,9)        ! B0(max_ht), b1..b8 per FVS species
    LOGICAL LGREG, GHAVE_DG(MAXSP), GHAVE_HG(MAXSP)
    REAL    GEMT, GTD, GELEV    ! per-STAND climate + elevation, set at stand setup

## GREGLOAD (pattern of GOMPLOAD)
- Activation: env FVS_GREGDG / FVS_GREGHG (or a GREGGROW keyword), with FVS_GREGDG_COEF / FVS_GREGHG_COEF
  pointing at config/greg_dg_coefficients.csv and greg_hg_coefficients.csv (header + SPCD n B0..). FIAJSP maps
  FIA SPCD to FVS species, as GOMPLOAD does.
- Per-stand climate: GEMT, GTD from a stand lookup (FVS_GREG_CLIMATE=path), keyed by stand id or fuzzed
  lat/lon, sourced from the projector's greg_emt_td_lookup.rds (export to CSV). GELEV from the stand record
  (PLOT/CONTRL ELEV, in feet, as Greg's forms expect). Resolve once per stand.

## GREGDG / GREGHG evaluators (validated; copy from test_gregdghg.f90)
    SUBROUTINE GREGDG(ISPC, DBH, CR, HT, BAL, ELEV, EMT, G)
      z = GDG(ISPC,1) + GDG(ISPC,2)*log((DBH+1)**2/(CR*HT+1)**GDG(ISPC,4))
            + GDG(ISPC,3)*BAL**GDG(ISPC,5)/log(DBH+2.7) + GDG(ISPC,6)*ELEV + GDG(ISPC,7)*EMT
      z = min(max(z,-30.),5.);  G = max(exp(z),0.)
    SUBROUTINE GREGHG(ISPC, HT, CR, CCFL, CCH, ELEV, TD, EMT, DHT)
      crp=max(CR,1e-4); cchp=max(CCH,0.); tdp=max(TD,0.)
      arg = -GHG(ISPC,2)*HT - GHG(ISPC,5)*CCFL - GHG(ISPC,9)*cchp**0.5 - GHG(ISPC,6)*ELEV
              + GHG(ISPC,7)*sqrt(tdp) + GHG(ISPC,8)*EMT
      DHT = GHG(ISPC,1)*GHG(ISPC,2)*GHG(ISPC,3)*crp**GHG(ISPC,4)*exp(arg)*(1.-exp(-GHG(ISPC,2)*HT))**(GHG(ISPC,3)-1.)
      DHT = max(DHT,0.)

## Tree-state inputs (availability)
- dbh: DBH array (have). cr: crown ratio ICR/100 (have). ht: HT array (have).
- bal: basal area in larger trees, BALMOD (have). ccfl: crown competition factor in larger (FVS computes).
- cch: crown closure at tip, already recomputed by GOMPCCH each cycle; reuse it (verify HG wants the 0-1
  fraction, not the gompit affine scale; pass the pre-affine ORGANON cch_hat fraction to HG).
- elev/EMT/TD: per stand (GELEV/GEMT/GTD from GREGLOAD).

## Validation plan (mirrors the mortality path)
1. Unit (done): standalone GREGDG/GREGHG reproduce the projector to 8e-8 / 4e-6; increments plausible.
2. Stand A/B: once wired + rebuilt, run net01.key baseline vs FVS_GREGDG/FVS_GREGHG and confirm dbh/ht
   trajectories track the projector NE greg metrics (conus_eq_ne_greg_metrics.csv).

## Status
Data + evaluators validated and committed (branch feat/greg-dg-hg-native-data). Remaining: author GREGMC.f90 +
GREGLOAD + the DGF/HTGF one-line hooks (maintainer-domain engine edit + rebuild) and the per-stand EMT/TD
lookup export. The cr update is the fvs-conus kernel (non-Greg), so option 1 leaves CR to the native/our kernel.
