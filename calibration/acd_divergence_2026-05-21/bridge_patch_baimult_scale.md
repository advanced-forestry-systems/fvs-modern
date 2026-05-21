# Bridge patch: decouple Acadian calibration from FVS baimult, fix the scale

This is a ready to apply patch for the two calibration issues in Driver A. It is
written as before and after snippets rather than a blind in place edit, because
the full pipeline (fvsOL plus the variant DLL) must be exercised before the
change ships. Verify with `bridge_patch_verify.R` (the scale math) and
`acd_layer2_attribution.R` (the model behaviour).

A ready to apply form of Change 1 is provided as `make_fvs_calib_acd.patch`,
verified to apply cleanly (`patch -p1 --dry-run`), preserve the file's CRLF line
endings, and parse. From the Interface checkout root:

    patch -p1 < calibration/acd_divergence_2026-05-21/make_fvs_calib_acd.patch

It inserts an override block at the end of `make_fvs_calib()` that reads the
calibration csv (path from `FVS_ACD_CALIB` or the packaged `fvsOL/extdata` copy)
and replaces `dDBH.mult`, `dHt.mult`, `mort.mult`, falling back to the FVS values
if the csv is absent. Note: the fitted `dDBH.mult` is already a diameter scale
multiplier, so the line 2232 application stays as is; the basal area conversion
in Change 2 is only needed if you instead honor a raw FVS `baimult`.

## Change 1 (main): source the Acadian calibration table, not FVS baimult

The Acadian model should carry its own annual per species calibration, fit by
`annualized_calibration.R`, instead of borrowing the FVS-NE `baimult`,
`htgmult` and `mortmult`. Produce `acd_annual_calibration.csv` with columns
`SP, dDBH.mult, dHt.mult, mort.mult`, then in `make_fvs_calib()`
(AcadianGY.R near line 1327) read that table instead of the FVS attributes.

Before (AcadianGY.R, inside make_fvs_calib):

    calib.fvs = fvsGetSpeciesAttrs(c("baimult","htgmult","mortmult", ...))
    ... rename(dDBH.mult = baimult, dHt.mult = htgmult, mort.mult = mortmult) ...

After (read Acadian annual calibration, keep FVS only for size caps):

    acd.cal = read.csv(system.file("extdata","acd_annual_calibration.csv",
                                   package="fvsOL"), stringsAsFactors = FALSE)
    # keep maxdbh / maxht from FVS for the size caps, but take the growth and
    # mortality multipliers from the Acadian table
    calib.fvs = fvsGetSpeciesAttrs(c("maxdbh","maxht","minmort","maxdbhcd"))
    ... join acd.cal by SP, coalesce dDBH.mult / dHt.mult / mort.mult to 1 ...

This removes the conceptual mixing of two calibration systems and is the lever
that actually matters (per species effects up to about 9 percent in the
attribution run).

## Change 2 (secondary): if baimult is still honored, convert through basal area

If you keep an option to honor a user supplied `baimult` (a basal area increment
multiplier), do not apply it to the linear diameter increment. Convert it.

Before (AcadianGY.R, in AcadianGYOneStand, line 2232):

    dDBH = dDBH * dDBH.thin.mod * dDBH.SBW.mod * dDBH.form.risk.mod * dDBH.mult

After (apply the legitimate linear modifiers, convert the basal area multiplier):

    dDBH = dDBH * dDBH.thin.mod * dDBH.SBW.mod * dDBH.form.risk.mod
    # dDBH.mult holds a basal area increment multiplier; convert so the realized
    # basal area increment is exactly dDBH.mult times the unscaled increment
    dDBH = sqrt(DBH^2 + dDBH.mult * ((DBH + dDBH)^2 - DBH^2)) - DBH

The thinning, spruce budworm and form/risk modifiers are genuine diameter
increment ratios, so they stay as linear multipliers. Only `dDBH.mult` (the
basal area multiplier) needs the conversion. Note the small inside bark versus
outside bark caveat: FVS DDS is inside bark, the Acadian increment is outside
bark, so the conversion is exact for the outside bark basal area target and very
close for the inside bark one.

## Change 3 (hygiene): align option defaults and remove the dead mortType

In `customRun_fvsRunAcadian.R`, make the wrapper defaults match the standalone
defaults you validate against (INGROWTH, CutPoint), and remove the
`mortType="continuous"` ops entry (line 214) because `AcadianGYOneStand` forces
discrete in 12.3.5, so the entry is silently ignored and only invites drift.

## Change 4 (separate effort): CSI and HT/CR initialization

Decide CSI once (push the Acadian CSI into the Event Monitor so the wrapper
reads it rather than re-deriving from site index, customRun lines 169 to 171),
and re-enable the Acadian HT and CR dubbing (customRun lines 95 to 140) once the
stop point 7 issue is resolved, or guarantee every input tree carries measured
HT and CR. These are tracked in the diagnosis memo as Drivers B and C.
