# AcadianGY 12.3.8: ingrowth STAND/PLOT inheritance (#127, part 2 of 2)

2026-05-26. 12.3.8 = 12.3.7 plus a one-block fix to the ingrowth path. It is a
strict improvement: with INGROWTH off it is identical to 12.3.7/12.3.6/12.3.5.
Multi-cycle harnesses now retain recruits across cycles.

## The bug (fixed)

Part 1 (12.3.7) made recruits grow by giving them neutral calibration
multipliers and species size caps before `bind_rows(tree, ingrow)`. That was
necessary but not sufficient. Recruits still inherited `STAND=1, PLOT=1` from
`ING.TreeList`'s default Sum.temp row instead of the parent stand's STAND/PLOT.

Direct evidence on FIA plot 243883099489998 with 12.3.7, CutPoint=0, one cycle:

    unique STAND values in OUTPUT frame:
                  1 243883099489998
                 12              29

29 survivors kept the FIA PLT_CN. The 12 recruits all got STAND=1. On the next
cycle, any multi-stand harness that iterates over `unique(trees$STAND)` and
looks up stand-level inputs by STAND finds no match for "1", the per-stand
`AcadianGYOneStand` call errors with "missing value where TRUE/FALSE needed,"
and the recruits are silently dropped.

That is exactly why the original v21 result on 100 FIA plots over 10 years was
byte-identical between baseline (CutPoint 0.95, no ingrowth) and ingrowth_fix
(CutPoint 0): the 12.3.7 carry-through fix made recruits, but the harness
discarded them on every subsequent cycle. The earlier INGROWTH_FIX_12.3.7.md
hypothesis (a NA ingrowth-rate covariate on FIA stands) is FALSE. The
instrumented probe shows BAPH, tph, qmd, pHW.ba all have NA count = 0 and IPH
in 5..17 on every FIA stand.

## The fix

In `AcadianGYOneStand`, extending the existing 12.3.7 ingrow block (right
before `tree = dplyr::bind_rows(tree, ingrow)`):

    if (!is.null(ingrow) && nrow(ingrow) > 0) {
      .cap <- unique(tree[, c("SP","max.dbh","max.height")])
      .cap <- .cap[!duplicated(.cap$SP), ]
      .mi  <- match(ingrow$SP, .cap$SP)
      ingrow$dDBH.mult  <- 1
      ingrow$dHt.mult   <- 1
      ingrow$mort.mult  <- 1
      ingrow$max.dbh    <- ifelse(is.na(.cap$max.dbh[.mi]),    200, .cap$max.dbh[.mi])
      ingrow$max.height <- ifelse(is.na(.cap$max.height[.mi]),  60, .cap$max.height[.mi])
      # 12.3.8: STAND/PLOT inheritance
      if ("STAND" %in% names(tree) && length(unique(tree$STAND)) == 1) {
        ingrow$STAND <- unique(tree$STAND)
      }
      if ("PLOT" %in% names(tree) && length(unique(tree$PLOT)) == 1) {
        ingrow$PLOT <- unique(tree$PLOT)
      }
    }

`AcadianGYOneStand` operates on one stand per call, so survivors share a single
STAND value. Recruits now inherit that value and stay with the parent across
cycles. Version tag bumped to AcadianV12.3.8.

## Unit proof

Same plot, same cycle, same options, 12.3.8 instead of 12.3.7:

    unique STAND values in OUTPUT frame:
    243883099489998
                 41

All 41 trees (29 survivors + 12 recruits) carry one STAND value.

## Harness proof

10 FIA plots, 3 annual cycles, CutPoint=0, MORTCAL on, sourcing 12.3.8:

| config                          | TPA     | BA bias %    | QMD (cm) |
|---------------------------------|---------|--------------|----------|
| baseline (CutPoint 0.95)        | 1013.47 | +4.73        | 4.956    |
| ingrowth_fix (12.3.8, CP=0)     | 1025.81 | +4.66        | 4.910    |

The configurations now produce different aggregates. Same harness with 12.3.7
returned byte-identical numbers for both. No `[V23 err]` lines and no phantom
`sid=1` calls in the trace.

## Status

Parts 1 (frozen recruits, 12.3.7) and 2 (recruit STAND identity, 12.3.8) are
fixed. The full 100-plot 10-year FIA validation (v24) is the next confirmation
that the QMD residual closes from 5.15 toward observed 4.97 and BA bias drops
below +10.9%; that run is queued on Cardinal.
