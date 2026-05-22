# AcadianGY 12.3.6 validation (in-source opt-in MORTCAL)

2026-05-21. 12.3.6 is canonical 12.3.5 plus the size-dependent mortality
correction (task #126b) integrated natively into the per-cycle mortality step of
`AcadianGYOneStand`, replacing the earlier external `AcadianGY_12.3.5_mortcal.r`
wrapper.

## What changed from 12.3.5 (62 lines)

- Adds `.mortcal_surv_ratio(dbh_in)` (observed/AcadianGY 5-yr survival ratio by
  post-growth DBH class, Maine FIA v10 diagnostic).
- In the per-cycle tree update, when `ops$MORTCAL = TRUE`, multiplies SURVIVOR
  EXPF by `surv_ratio(DBH)^(1/MORTCAL_INTERVAL)` (default interval 5) BEFORE
  ingrowth is bound in.
- Default is off: with `ops$MORTCAL` not TRUE the model is identical to 12.3.5.

Improvement over the wrapper: the wrapper trimmed all trees post-step including
fresh ingrowth; 12.3.6 haircuts survivors only, so ingrowth is preserved.

## Validation (test_12_3_6.R, Cardinal, R 4.4.0)

Sourced canonical 12.3.5 and 12.3.6, projected the same 60-tree mixed Acadian
stand 25 years (INGROWTH off to isolate mortality):

| model | BA m2/ha | TPH |
|---|---|---|
| 12.3.5 canonical | 25.397 | 447.69 |
| 12.3.6 MORTCAL off | 25.397 | 447.69 |
| 12.3.6 MORTCAL=TRUE | 20.017 | 357.60 (BA -21.2%) |

A. Equivalence: 12.3.6 with MORTCAL off is identical to canonical 12.3.5 (BA and
   TPH match exactly). B. Direction: MORTCAL=TRUE removes survivors and lowers
   BA/TPH as designed. Both pass.

## Notes and one bug to fix

- Version tag bug: `AcadianVersionTag` is still `"AcadianV12.3.5"` in the 12.3.6
  file; it should read `"AcadianV12.3.6"` so runs are self-identifying.
- Condition-specific (file header, #128): the survival ratios are Maine FIA
  calibrated and OVER-correct on Canadian CFI (MAGPlot NB), where basal area is
  already unbiased. This matches our MAGPlot finding (R AcadianGY BA bias ~0% on
  262 NB plots). So keep MORTCAL OFF for the MAGPlot / Canada Fortran runs.
- Partial (#127): on 200 Maine FIA plots MORTCAL roughly halves the BA
  over-projection (+15.4% to +8.6%) but the residual is a QMD / small-tree
  ingrowth composition effect, not mortality. The next fix is in the
  recruitment/ingrowth submodel.

## Cross-validation: FIA vs Canadian CFI (the #128 caveat, empirical)

The in-source 12.3.6 was run on both datasets, off and on, completing the 2x2
matrix (FIA in-source = Aaron's v16; MAGPlot in-source = `cardinal_magplot_insource_v16.R`):

| dataset | MORTCAL off | MORTCAL on (in-source) |
|---|---|---|
| Maine FIA (200 plots) | BA +15.4%, R2 0.417 | BA +8.6%, R2 0.483 |
| Canadian NB MAGPlot (262 pairs) | BA -0.04%, R2 0.878 | BA -6.51%, R2 0.869 |

Reading: on Maine FIA the correction halves a real over-projection (+15.4 to
+8.6). On Canadian CFI, where AcadianGY is already unbiased (-0.04), the same
correction drives basal area to -6.51 (over-correction). This empirically
confirms #128: the correction is FIA-specific, not a universal model bias. The
MAGPlot off case (-0.04, R2 0.878) also confirms 12.3.6 with MORTCAL off
reproduces canonical 12.3.5 on real NB data. The in-source on (-6.51) is a hair
closer to zero than the wrapper's -6.56, consistent with ingrowth being
preserved.

## Recommendation

Adopt 12.3.6 as `fvsOL/inst/extdata/AcadianGY.R` in the Interface checkout
(replacing 12.3.5). The version tag has been corrected to `AcadianV12.3.6` in
`AcadianGY_12.3.6_tagfixed.r`. It is a safe drop-in (off by default). Enable
`ops$MORTCAL=TRUE` only for FIA-like Maine conditions, never for Canadian CFI /
MAGPlot (it over-corrects there, as the matrix shows).
