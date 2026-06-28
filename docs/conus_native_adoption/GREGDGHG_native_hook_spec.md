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

## Engine substitution points (NE; verified against the source -- read carefully, subtle)
Both substitutions are guarded per species by a GHAVE flag, like morts.f90 (IF(LGREG.AND.GHAVE_DG(ISPC))...),
so unfit species keep the native equations.

### DGF (diameter growth) -- the WK2 it stores is LOG of inside-bark DDS, not raw DDS
The native (dgf.f90) iterates an ANNUAL DG 10 times (DO 1000 ILOOP=1,10) updating TEMD(I), then computes:
    BARK  = BRATIO(ISPC, TEMD(I), HT(I))
    DIAGR = (TEMD(I) - DIAM(I)) * BARK            ! inside-bark period increment
    DDS   = DIAGR * (2.0*DIAM(I)*BARK + DIAGR)
    WK2(I)= ALOG(DDS) + COR(ISPC)
So the Greg hook must produce the SAME WK2 form. For a covered tree, replace the 10-iteration native block with
a 10-step Greg loop (cr, ht, bal, elev, emt held fixed across the cycle; only D updates):
    D = DIAM(I); do y=1,FINT_years:  D = D + GREGDG(ISPC,D,cr,ht,bal,elev,emt)   ! outside-bark dbh
    BARK  = BRATIO(ISPC, D, HT(I))
    DIAGR = (D - DIAM(I)) * BARK
    DDS   = DIAGR * (2.0*DIAM(I)*BARK + DIAGR);  if (DDS < 1e-6) DDS = 1e-6
    WK2(I)= ALOG(DDS)            ! OMIT COR(ISPC): that is the native NE-TWIGS calibration, not Greg's
Greg's dg is an OUTSIDE-bark dbh increment (FIA dbh is outside bark); BRATIO converts to the inside-bark DDS
the engine consumes. Inserts just before the WK2(I)=ALOG(DDS)+COR(ISPC) line, replacing it for covered species.

### HTGF (height growth) -- bypass the native modifiers for Greg
The native (htgf.f90) builds HTG(I) from HTCALC then applies BALMOD, a relative-height temper, OLDRN random,
SCALE, XHT and HTCON. Greg's HG already encodes competition (ccfl, cch) and climate, so for a covered tree
REPLACE the whole native HTG(I) with:
    HTG(I) = GREGHG(ISPC, HT(I), ICR(I)/100., ccfl, cch_frac, elev, td, emt) * FINT_years
Insert after the native HTG(I) is first set, guarded by GHAVE_HG; do NOT then apply BALMOD/SCALE/HTCON.

### Per-tree inputs and the two gotchas
- dbh = DIAM(I); cr = ICR(I)/100.0; ht = HT(I) (ARRAYS common). FINT_years from CONTRL.
- bal (ft2/ac, basal area in larger): computed by BADIST/BALMOD; confirm the per-tree BAL array name/common
  (BADIST stores it) -- needed raw, not the BAGMOD modifier.
- ccfl (crown competition factor in larger): FVS computes a CCF-in-larger; confirm the per-tree array.
- GOTCHA 1 (cch scale): GREGHG wants the ORGANON crown-closure FRACTION (0-1). GOMPCCH currently fills CCHT(I)
  on the gompit AFFINE scale (0.062 + 0.0036*cch_hat). Expose the pre-affine cch_hat (0-1) for HG, or recompute;
  do NOT feed the affine CCHT into GREGHG.
- GOTCHA 2 (bark): DG must round-trip through BRATIO as above; feeding Greg's outside-bark increment straight
  into DDS without BARK would be wrong.
- elev/emt/td: per stand. FIRST VERSION supply as scalars via env (FVS_GREG_EMT, FVS_GREG_TD) + stand ELEV, so
  the single-stand net01.key A/B is runnable; FULL version reads the per-stand greg_emt_td_lookup (export to CSV).

These mechanics (LOG-DDS + BRATIO, modifier bypass, cch fraction, BAL/CCFL arrays) are why this is a careful
maintainer-domain edit with a rebuild-and-A/B loop, not a blind patch. The evaluators are already validated;
the risk is entirely in the wiring above.

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
