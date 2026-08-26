# GREGHTDBH: native engine hook for the Marshall CONUS height-diameter model

Design and verification record for wiring the Marshall CONUS height-diameter (HT-DBH) fits into the FVSne
growth loop, mirroring the working Greg DG hook (base/gregdghg.f90 + GREGMC common + env-var activation + a
guarded one-line substitution in the host subroutine). This hook substitutes the predicted total height from
DBH used to dub missing tree heights, replacing the native Wykoff / Curtis-Arney prediction for covered
species while leaving uncovered species and all diameter-from-height calls unchanged.

## Where it substitutes

FVSne dubs missing heights from diameter in `ne/htdbh.f90`, called from `ne/cratet.f90` in MODE 0 (diameter
provided, height desired). The native routine uses either a Wykoff form (`IWYKCA==0`) or a Curtis-Arney form
(`IWYKCA==1`), keyed by the FVS species index. The hook is inserted at the end of the MODE 0 branch, after the
native `H` is computed:

    IF (LGREGHTDBH) THEN
      IF (GHAVE_HD(ISPC)) THEN
        GHDVAL = GREGHTDBH(ISPC, D)
        IF (GHDVAL .GT. 0.0) H = GHDVAL
      ENDIF
    ENDIF

MODE 1 (height provided, diameter desired; used by REGENT) is deliberately left on the native inverse. The
substitution is centralized in HTDBH, so both CRATET call sites (lines 353 and 425) are covered by one edit.

`GREGLOADHD` is called once per stand from CRATET (right after the CRATET DBCHK, before any HTDBH dubbing) so
that `LGREGHTDBH` and the per-species coverage flags are always defined before HTDBH reads them. The loader
carries a `SAVE :: LDONE` guard, so the file is read only once per run.

## Activation (no recompile to toggle)

Two environment variables, read once in `GREGLOADHD`:

    FVS_GREGHTDBH        1/on/true to enable the CONUS HT-DBH substitution. Default OFF and silent when unset
                         (blank, 0, n, or N leave LGREGHTDBH=.FALSE. and return immediately).
    FVS_GREGHTDBH_COEF   path to conus_htdbh_coefficients.csv. If empty or unopenable, the hook logs a one-line
                         message to JOSTND and stays disabled.

When disabled the engine behavior is bit-for-bit the native FVSne. When enabled, only species carrying a
fitted row are affected; every other species keeps the native Wykoff / Curtis-Arney height.

## The six forms

Transcribed verbatim from the `model.formula` strings in the Marshall
`sf_integration_dev/htdbh_native/marshall_fits/CONUS_HDfit_Model_*.CSV` files. `bh = 4.5` ft; DBH in inches;
HT returned in feet. `GREGHTDBH` dispatches on the per-species `model_id` (1..6):

    model 1: HT = bh + exp(B1 + B2*DBH^(-1.0))
    model 2: HT = bh + exp(B1 + B2*DBH^B3)
    model 3: HT = bh + B1*(1.0 - exp(B2*DBH))^B3        (Chapman-Richards)
    model 4: HT = bh + exp(B1 + B2/(DBH + 1.0))
    model 5: HT = bh + exp(B1 + B2/(DBH + B3))
    model 6: HT = bh + B1*(1.0 - exp(B2*DBH^B3))        (Weibull)

B3 is unused for models 1 and 4 (blank in the CSV, parsed to 0.0). `GREGHTDBH` returns a sentinel of -1.0 for
uncovered species, for non-positive DBH, and for an out-of-range model_id, so the caller keeps the native
value in every degenerate case.

## Coefficient schema

`conus_htdbh_coefficients.csv` (staged on the branch at sf_integration_dev/htdbh_native/):

    SPCD,model_id,B1,B2,B3,nobs,AIC,RMSE,meanBIAS,ierr,tier

The loader reads only the first five columns (SPCD, model_id, B1, B2, B3). SPCD is the FIA species code; it is
matched to the FVS species index via FIAJSP (as GREGLOADDG does). A dedicated field-tolerant parser
(GREGHD_PARSE5) splits on commas and treats a blank B3 field as 0.0, avoiding list-directed READ swallowing
across the empty field. The staged table covers 412 species (392 tier 1, 20 tier 2); the remaining CONUS
species have no row and fall back to the native equations.

## New / changed files

    src-converted/common/GREGHDMC.f90    new common block: LGREGHTDBH, NGREGHD, GHDMOD(MAXSP),
                                         GHD(MAXSP,3), GHAVE_HD(MAXSP)
    src-converted/base/greghtdbh.f90     new: GREGLOADHD loader, GREGHD_PARSE5 parser, REAL FUNCTION GREGHTDBH
    src-converted/ne/htdbh.f90           wired: GREGHDMC include, GREGHTDBH declaration, guarded MODE 0
                                         substitution
    src-converted/ne/cratet.f90          wired: CALL GREGLOADHD once per stand before height dubbing
    src-converted/bin/FVSne_sourceList.txt  added ../base/greghtdbh.f90 after gregdghg.f90

## Verification

Standalone harness `sf_integration_dev/htdbh_native/test_greghtdbh.f90` transcribes the same six forms and,
compiled real(kind=8) (double precision; float32 only reaches ~4e-6, cf. PR #95
fix/gregdghg-real8-precision), evaluates every covered species over a DBH grid of 1..50 in, writing
`htdbh_fortran.csv`. The R script `validate_greghtdbh.R` re-evaluates each species' exact model.formula into
`htdbh_reference_R.csv` and compares. Result (412 species x 50 DBH = 20,600 rows):

    model 1 : n= 1600  max_abs_diff = 4.974e-13
    model 2 : n= 4050  max_abs_diff = 5.116e-13
    model 3 : n= 2950  max_abs_diff = 5.116e-13
    model 4 : n= 2800  max_abs_diff = 4.974e-13
    model 5 : n= 6700  max_abs_diff = 4.974e-13
    model 6 : n= 2500  max_abs_diff = 4.974e-13
    OVERALL max abs diff = 5.116e-13   (PASS, < 1e-6)

The full FVSne executable builds cleanly with the hook (deployment/scripts/build_fvs_executables.sh ne: 557
objects, 1 built, 0 failed, 7.9M FVSne), and the symbols gregloadhd_, greghtdbh_, and greghd_parse5_ are
present in the linked binary.
