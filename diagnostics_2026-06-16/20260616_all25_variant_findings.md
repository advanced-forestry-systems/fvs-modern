# All-variant stress test: does calibration improve the FVS engine? (2026-06-16)

21 of 25 variants run through the FVS engine on FIA remeasurement (AK, and the Canada variants BC/ON,
have no usable FIA here). Default vs calibrated, BA and density bias vs observed. Figure:
`all25_engine_2026-06-16.png`; table: `all25_engine_default_vs_calibrated.csv`.

## Headline results

1. FVS over-predicts basal area in 20 of 21 variants (the lone exception is EC, which under-predicts
   at -7%). The over-prediction ranges from +5% (WS) to +42% (EM). This matches the well-known FVS
   behavior, and it confirms the engine benchmark is faithful (unlike the R-equation harness, which
   under-predicts because it omits ingrowth and stand dynamics).

2. The calibrated config equals default within about one point in 16 of 21 variants: the per-species
   multipliers are inert at the stand level (they are centered on 1, so they redistribute among
   species and net to zero).

3. Three variants (IE, CI, KT) show a lower calibrated BA bias (e.g. IE 19.8 to 8.4). But their
   density crashes at the same time (IE TPH -24.9 to -36.7, CI -13.8 to -26.6, KT -28.9 to -38.7).
   The BA "improvement" is OVER-THINNING, the same compensating artifact as the original SDIMAX bug:
   killing trees lowers BA toward observed. It is driven by large mortality multipliers in those
   variants, not by better growth. So it is not a legitimate improvement.

Net: calibration does not legitimately improve the FVS engine's basal-area prediction in ANY of the 21
variants. Wherever it appears to help, it is over-thinning.

## Why (answering "I do not understand how calibration would not improve predictions")

The calibration improved the EQUATIONS, and applied directly (the R `fia_benchmark`) it does reduce
bias. The problem is two-fold and both are in delivery, not in the fitting:

- The multipliers are RELATIVE. `dds_multiplier = exp(b0_species - mean_b0)`, centered on 1. The
  differencing against the mean discards the absolute level, which is exactly the term that would
  correct an over-prediction. So they reshuffle species, not the stand-level bias.
- The refit EQUATIONS never reach the engine. The engine runs native FVS plus these inert multipliers;
  the actual refitted Wykoff/ORGANON/Chapman-Richards equations are only injected in the prototype that
  does not load trees (the fvs2py blocker).

And the validation that showed improvement (the R harness) UNDER-predicts, the opposite of the real
engine, so it could not have surfaced this: tuning toward zero from below is the wrong direction for an
engine that sits above.

## What would actually improve it

A SIGNED correction that pulls growth down where FVS over-predicts (most variants) and up where it
under-predicts (EC, and the Western Sierra per Batista et al.). Demonstrated on NE: a global
diameter-growth slowdown (BAIMULT < 1) reduces BA and QMD bias monotonically. The principled version is
to inject the refitted equations (carrying their absolute level), once the fvs2py tree-loading is
fixed. Either way the target is the absolute stand-level bias, not per-species relative ratios.

## On uncertainty (status)

The predictive-interval fix is VERIFIED (parameter-only intervals cover ~13%; adding residual variance
restores ~93%), but it is not yet deployed in the production `21_uncertainty_propagation.R`. So
uncertainty is now correctly understood and the fix is proven, but "properly accounted for" requires
the one bounded code change.
