# Test of the updated FVS-ACD R code (AcadianGY_12.3.5_mortcal.r)

2026-05-21. The updated code adds an opt-in size-dependent mortality correction
to the canonical AcadianGY 12.3.5. This is the fix for the issue flagged in the
diagnosis: `mort.mult` is plumbed but never applied in the canonical model, so
its mortality levers cannot impose the observed size-by-class survival. The
correction instead applies, each annual cycle, a post-step EXPF multiplier
`ratio(DBH class)^(1/interval)`, where `ratio` is the observed/AcadianGY 5 year
survival ratio by DBH class from the FIA diagnostic.

## What was tested

`test_mortcal.R` sources the updated file (which sources the canonical model and
defines `AcadianGYOneStandMortCal`) and projects a 60 tree mixed Acadian stand
annually, comparing the canonical model, the wrapper with `MORTCAL=FALSE`, and
the wrapper with `MORTCAL=TRUE`. The uploaded file is byte identical to the copy
on Cardinal (`seven_islands/AcadianGY_12.3.5_mortcal.r`).

## Results

| horizon | config | BA m2/ha | TPH |
|---|---|---|---|
| 5 yr | canonical | 17.14 | 447.7 |
| 5 yr | MORTCAL=FALSE | 17.14 | 447.7 |
| 5 yr | MORTCAL=TRUE | 16.36 | 427.7 (BA -4.5%) |
| 25 yr | canonical | 25.40 | 447.7 |
| 25 yr | MORTCAL=FALSE | 25.40 | 447.7 |
| 25 yr | MORTCAL=TRUE | 20.02 | 357.6 (BA -21.2%) |

Both design guarantees hold:

1. Reversibility. With `MORTCAL` not TRUE the wrapper is identical to the
   canonical model (BA and TPH match exactly at both horizons).
2. Direction. With `MORTCAL=TRUE` the correction removes stems (all survival
   ratios are below 1), lowering TPH and BA, which is the intended removal of the
   basal area over-projection. The effect compounds with horizon because the
   annual ratio is applied every year (over 25 years EXPF is scaled by roughly
   the 5 year ratio to the fifth power).

This confirms the updated R code behaves as specified. The remaining residual
(about +8.6 percent on the 200 plot Maine FIA set per the file header) is a
diameter growth and ingrowth composition issue, not mortality, consistent with
the diagnosis.

## How it composes with the annual diameter/height calibration

The mortcal correction is orthogonal to `acd_annual_calibration.csv`: the
calibration table scales diameter (and provisionally height) increments through
`dDBH.mult`/`dHt.mult`, while mortcal scales survival through EXPF. They can be
applied together by running `AcadianGYOneStandMortCal` with `ops$MORTCAL=TRUE`
and the calibration table loaded via the make_fvs_calib patch.
