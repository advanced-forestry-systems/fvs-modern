# Four-arm comparison: default / keyword-calibrated / fvs-conus / combined
2026-06-18. Item 2 of the integration roadmap. Establishes the shared headline comparison and resolves the
projector-vs-engine sign discrepancy by reporting each approach as its within-framework improvement.

## The framework problem and how it is resolved

The four arms do not all live in one projection machinery yet:

- Arms A (default) and B (keyword-calibrated: brms max SDI + density-dependent recruitment + BAIMULT) run
  in the FVS engine, where the keyword levers apply. The engine over-predicts undisturbed BA (about +9%).
- Arm C (fvs-conus species-free equations) runs in the fvs-conus standalone projector, which under-predicts
  undisturbed BA (about -12%). Opposite sign, different machinery.
- Arm D (combined: fvs-conus equations + brms max SDI + density-dependent recruitment) needs the fvs-conus
  equations running inside the FVS engine. That is the fvs2py in-process tree-loading injection, the one
  documented maintainer-level blocker.

A single absolute-bias four-arm therefore conflates the equation swap with the framework offset. The fix
used here is to report each arm as the improvement relative to its own framework's default. The
within-framework |bias| reduction is framework-invariant and is the honest comparable quantity until the
fvs2py injection lands and lets all four arms share the engine.

## NE, disturbance-clean basis (both frameworks)

Engine arms A and B are out-of-sample on spatial fold B (county-hash folds, calibration derived on fold A
only), n = 128, all four metrics with bootstrap CIs. Projector arms (default and C) are the 21,811
COND-undisturbed NE conditions from the 2026-06-17 disturbance-clean run. The all-variant engine run
(ne, sn, kt, pn, cr, ut, nc, ec, wc; job 11745221) is completing and will populate a cross-variant median
table; NE is shown here and is complete in both frameworks.

| metric | A default (engine) | B keyword-cal (engine) | engine improvement | default (projector) | C fvs-conus (projector) | projector improvement |
|---|---|---|---|---|---|---|
| BA  | +8.7 | +8.5 | 0.2 | -12.3 | -7.2 | 5.1 |
| TPH | -14.6 | -11.2 | 3.4 | -7.2 | -6.2 | 1.0 |
| QMD | +12.1 | +2.1 | 10.0 | -5.3 | -3.1 | 2.2 |
| VOL | +21.2 | +21.2 | 0.0 | -15.7 | -9.2 | 6.5 |

(Engine values are % bias, default to calibrated; engine CIs in fourarm_engine_20260618.csv.)

## The finding: the two efforts are complementary, not competing

The within-framework improvements separate cleanly by metric:

- Keyword calibration (engine, arm B) delivers the size and density gains: QMD bias collapses (+12.1 to
  +2.1) and TPH improves (-14.6 to -11.2). It does almost nothing for stand-level BA or merch volume,
  because brms max SDI and recruitment do not change per-tree growth level.
- fvs-conus equations (projector, arm C) deliver the level and scatter gains: BA improves (12.3 to 7.2)
  and volume improves most of all (15.7 to 9.2), with lower RMSE and better species consistency. They do
  little for QMD or TPH because the equations do not touch stand density or recruitment.

This is the arm-D hypothesis stated precisely: combining the fvs-conus growth equations (which fix BA and
volume level) with the brms max SDI and density-dependent recruitment (which fix QMD and TPH) should
capture both columns of improvement at once. It is the single most valuable remaining test and the shared
headline figure for both manuscripts.

## Status of each arm

| arm | definition | framework | state |
|---|---|---|---|
| A | default FVS | engine | NE done OOS with CIs; all-variant run completing (job 11745221) |
| B | brms maxSDI + density-dependent recruitment + BAIMULT | engine | NE done OOS with CIs; all-variant run completing |
| C | fvs-conus species-free equations | projector | done on NE disturbance-clean; extend beyond NE |
| D | fvs-conus equations + maxSDI + recruitment | engine (needs fvs2py) | BLOCKED on in-engine fvs-conus injection |

## To produce the fully unified single-framework four-arm

Two routes to put C and D in the engine and retire the cross-framework caveat:

1. fvs2py in-process tree loading: load the FIA tree list into the engine in memory, swap the diameter
   growth, height-diameter, and crown calls for the fvs-conus equation evaluations, project, and read the
   same summary. This is the clean path and unblocks both C and D directly.
2. Per-species multiplier emulation: express the fvs-conus per-species diameter-growth and height-diameter
   predictions as engine BAIMULT and HT-DBH keyword multipliers (ratio of fvs-conus to default per species),
   so arm C becomes an engine keyword run and arm D = C + the density levers. Approximate but unblocks the
   comparison now; requires the per-species fvs-conus prediction table (banked posteriors) across variants.

Engine harness: `fourarm_engine.py`. Engine results: `fourarm_engine_20260618.csv`. Projector arm C:
`fourarm_clean.py` and the fvs-conus disturbance-clean run.
