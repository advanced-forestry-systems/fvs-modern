# AcadianGY 12.3.7: ingrowth carry-through fix (#127, part 1 of 2)

2026-05-22. 12.3.7 = 12.3.6 plus a one-block fix to the ingrowth path. It is a
strict improvement: with INGROWTH off it is identical to 12.3.6/12.3.5.

## The bug (fixed)

`ING.TreeList` builds recruited trees with only STAND/PLOT/YEAR/TREE/SP/DBH/HT/
HCB/EXPF/pHT/pHCB/Form/Risk. It does NOT set `dDBH.mult`, `dHt.mult`,
`mort.mult`, `max.dbh`, or `max.height`. After `bind_rows(tree, ingrow)` those
columns are NA on every recruit, so on the next cycle `dDBH = dDBH * dDBH.mult`
is NA, then `DBH = DBH + coalesce(dDBH, 0)` leaves DBH unchanged. Recruits were
added but FROZEN at the 3 cm recruitment diameter and never grew into the stand.

Unit proof (single 60-tree stand, CutPoint=0): under 12.3.6 the 11 recruits came
out with NA count = 11 for all five columns; under 12.3.7 the NA count is 0, so
`dDBH.mult = 1` and the recruits grow.

## The fix

In `AcadianGYOneStand`, immediately before `tree = dplyr::bind_rows(tree, ingrow)`:

    if (!is.null(ingrow) && nrow(ingrow) > 0) {
      .cap <- unique(tree[, c("SP","max.dbh","max.height")])
      .cap <- .cap[!duplicated(.cap$SP), ]
      .mi  <- match(ingrow$SP, .cap$SP)
      ingrow$dDBH.mult  <- 1
      ingrow$dHt.mult   <- 1
      ingrow$mort.mult  <- 1
      ingrow$max.dbh    <- ifelse(is.na(.cap$max.dbh[.mi]),    200, .cap$max.dbh[.mi])
      ingrow$max.height <- ifelse(is.na(.cap$max.height[.mi]),  60, .cap$max.height[.mi])
    }

Recruits inherit species size caps from the survivors and neutral multipliers.
Version tag bumped to AcadianV12.3.7.

## A second issue remains (part 2 of 2)

On the FIA validation harness, 12.3.7 with CutPoint=0 still shows no net ingrowth
(baseline and ingrowth_fix both +10.9%, TPA 1019, QMD 5.15), even though the
synthetic debug stand recruits 11 trees. So a separate, data-level suppression
exists in the FIA path: the most likely cause is an NA covariate (e.g. `pHW.ba`
or `qmd`) reaching `Ingrowth.FUN`, making `IPH = NA`; then
`TreeCon = Sum.temp[Sum.temp$IPH > 0, ]` selects zero rows (NA > 0 is NA), so no
recruits are produced and no error is raised. This is the next step: guard the
ingrowth-rate covariates against NA in `AcadianGYOneStand`'s ingrowth block (and
confirm the FIA stand summaries do not produce NA `pHW.ba`/`qmd`). Once recruits
appear on FIA stands, re-check whether QMD drops from 5.15 toward observed 4.97
and the BA residual shrinks.

## Status

Part 1 (frozen recruits) is fixed and unit-validated in 12.3.7. Part 2 (FIA
stands compute NA ingrowth rate) is localized with a clear hypothesis and is the
remaining work to actually move the stand-level QMD/BA.
