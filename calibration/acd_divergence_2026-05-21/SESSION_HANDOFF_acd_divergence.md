# Session handoff: AcadianGY.R vs FVS-ACD divergence + annualized calibration

Date: 2026-05-21
Repo: holoros/fvs-modern
Branch: acd-bridge-followup-2026-05-20 (pushed; HEAD 8cb4ac0)
Worktree on Cardinal: /users/PUOM0008/crsfaaron/fvs-modern-acdbridge
Deliverables: calibration/acd_divergence_2026-05-21/ (also mirrored in the
selected fvs-online folder)

## Question this session answered

Aaron uploaded AcadianGY_12.3.5.r and reported that the R model output does not
match FVS-ACD, having ruled out annual versus periodic stepping and ingrowth as
the primary causes. He asked to understand why, and to develop calibration
factors for the R model on an annual rather than periodic basis.

## What was found

1. The equations are not the cause. The fork's
   fvsOL/inst/extdata/AcadianGY.R is byte for byte identical to the uploaded
   12.3.5 reference (0 differing lines). The stand basal area constraint is
   present in both.

2. The gap is created in the bridge, customRun_fvsRunAcadian.R. Ranked drivers,
   with the diagnosis corrected after quantification:
   - A. FVS calibration multipliers injected (baimult, htgmult, mortmult) where
     standalone uses 1. Largest part is simply that they are non unity and come
     from the wrong model (FVS-NE), not the Acadian model. The additional scale
     error (baimult is a basal area multiplier applied to a linear diameter
     increment) is real but only about 1.5 percent at the annual step.
   - B. Height and crown ratio dubbing disabled in the bridge, so FVS dubs with
     NE logic while standalone uses Acadian HTPred and HCBPred.
   - C. Climate site index sourced differently (bridge maps FVS site index;
     standalone uses CSI or 12).
   - D and E. Option defaults (INGROWTH off, CutPoint 0.5, ignored continuous
     mortType) and once per cycle ingrowth versus annual.
   - Annual versus periodic confirmed NOT a driver: the bridge already loops
     annually inside each FVS cycle (cyclen = 1). Matches Aaron's read.

3. If the comparison target is the compiled FVSacd binary rather than the
   customRun path, that binary runs NE Fortran equations in this fork, so it is
   a different model entirely.

## What was built and verified

- annualized_calibration.R: closed form basal area to diameter conversion plus a
  per species annual multiplier optimizer. Self test recovers known factors to 5
  decimals. Runs in base R.
- acd_divergence_decomposition.R: base R growth core (exact Kuehne dDBH form and
  stand basal area constraint) for offline attribution.
- acd_layer2_attribution.R: the real AcadianGYOneStand run on Cardinal under
  standalone versus bridge configs. It reproduces the base R core to within
  rounding (25 yr: multipliers +2.9 percent total, red spruce +9.0, yellow birch
  -4.7; stand total hides species divergence).
- bridge_patch_baimult_scale.md and bridge_patch_verify.R: ready to apply patch
  (read the Acadian annual table instead of FVS baimult; convert baimult through
  basal area if honored; align options; CSI and dubbing notes) with a base R
  proof that the legacy linear application is off by at most about 1.5 percent.
- fit_acd_annual_calibration_v3.R and acd_annual_calibration.csv: ANNUAL per
  species diameter multipliers fit against 575,461 real ACD FIA remeasurement
  records (17,326 plots) by iterating the real Acadian dDBH_fun. v3 corrects a
  units bug (recorded BAL is ft2/acre but dDBH_fun expects m2/ha; the first pass
  fed it in about 4.35x too high) and splits competition into BAL.SW and BAL.HW.
  34 species fitted; diameter growth RMSE 1.081 to 1.008 cm (6.8 percent
  reduction). With the corrected scale most multipliers are modestly above 1
  (balsam fir 1.14, red spruce 1.12, white pine 1.10) with larger corrections
  for fast hardwoods (red oak 1.92, red maple 1.50, aspen 1.48).
- fit_acd_height.R: dHt.mult now fitted from measured HT_t1/HT_t2 by a ratio of
  trimmed means (robust to the 23 percent noisy negative growth records). Range
  0.30 to 2.28, median 1.14; provisional given height data quality. mort.mult
  stays 1.0 (Acadian mortality is stand level; survival fit is a separate task).
- make_acd_calib_from_table.R: tested base R helper that builds the bridge
  calib.spp table from the CSV (keeps FVS size caps), so the patch is concrete.
- acd_annual_calibration_NOTES.md: method, the units correction, assumptions.

## Git state

Four commits on acd-bridge-followup-2026-05-20, all pushed:
- e105adb diagnosis, calibration framework, decomposition, real attribution
- 76ab3b5 bridge patch + verification + corrected memo
- 8cb4ac0 fitted annual calibration table + fit script + notes
(plus the prior 1f0068f from the earlier ACD calibrated row work)

## Remaining steps (for Aaron or next session)

1. Open the follow-up PR (needs your authenticated GitHub; gh on Cardinal is not
   logged in this session):
   https://github.com/holoros/fvs-modern/pull/new/acd-bridge-followup-2026-05-20

2. Apply the bridge patch in customRun_fvsRunAcadian.R / make_fvs_calib so the
   model reads acd_annual_calibration.csv instead of FVS baimult, then run the
   real customRun on a few stands to confirm R and FVS-ACD now agree.

3. Fit mort.mult on STATUSCD survival (the Acadian mortality is stand level, so
   this needs a survival model rather than the diameter machinery). Diameter and
   height multipliers are both in the table now; height is provisional given the
   noisy measured height data.

4. Decide CSI once (push Acadian CSI into the Event Monitor) and re-enable the
   Acadian HT and CR dubbing once the stop point 7 issue is resolved.

## How to reproduce on Cardinal

    ssh crsfaaron@cardinal.osc.edu        # OSC key id_osc
    cd /users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/acd_divergence_2026-05-21
    module load gcc/12.3.0 R/4.4.0
    Rscript annualized_calibration.R          # self test
    Rscript acd_divergence_decomposition.R    # offline attribution
    Rscript acd_layer2_attribution.R          # real model attribution
    FIT_CAP=6000 Rscript fit_acd_annual_calibration.R   # refit the table

The unrelated debug task from the earlier thread (task 131, the ACD calibrated
row in the FIA benchmark) is still open and independent of this work.
