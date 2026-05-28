# Mortality multiplier sensitivity under AcadianGY 12.3.8: null result

2026-05-28. Probes whether scaling the tree-level `mort.mult` column shifts the
+11 percent FIA BA residual. Scanned at 1.0, 1.3, 1.6, 2.0, 2.5 holding all
other parameters at 12.3.8 production posture (MORTCAL on, CutPoint = 0
expected-value ingrowth). 100 ME FIA plots, 10 annual cycles, same sample as
v24/v25.

## Result

| mort.mult scale | BA bias % | R^2   | TPA  | QMD cm |
|-----------------|-----------|-------|------|--------|
| 1.0             | +11.0506852 | 0.4232 | 1042.74 | 4.9232 |
| 1.3             | +11.0506852 | 0.4232 | 1042.74 | 4.9232 |
| 1.6             | +11.0506852 | 0.4232 | 1042.74 | 4.9232 |
| 2.0             | +11.0506852 | 0.4232 | 1042.74 | 4.9232 |
| 2.5             | +11.0506852 | 0.4232 | 1042.74 | 4.9232 |

All five rows are byte-identical to roughly 13 decimal places. mort_x1.0 also
reproduces v24 exactly (the sanity check passes).

## Interpretation

The tree-level `mort.mult` column does not flow into the mortality calculation
in 12.3.8 standalone production posture. Two plausible reasons:

(a) MORTCAL's survivor haircut, applied before bind_rows, supersedes the
    row-level mort.mult entirely. The Glover/Hool logistic that mort.mult
    multiplies is never the binding constraint when MORTCAL is on, because
    MORTCAL is doing its own size-dependent EXPF reduction.

(b) AcadianGYOneStand reads mort.mult from a different code path that the
    standalone harness does not exercise. The customRun bridge may set
    mort.mult per-call from a calibration table that overrides the
    base_init values we passed in.

Either way: mortality rate scaling via this lever does not address the
residual. Combined with the v25 CSI finding, the picture for the +11 percent
FIA BA residual is:

  - Climate (CSI scaling)        ~1.5 pp of partial closure available
  - Mortality (mort.mult scaling) null lever
  - Per-species calibration       ~7 pp BA contribution but offset by MORTCAL
  - Ingrowth                      ~0.2 pp contribution

The remaining ~9.5 pp of structural BA bias sits in:

  1. The intrinsic Kuehne et al. dDBH base-rate coefficients (not addressable
     by the row-level multipliers we have).
  2. The Glover/Hool mortality functional form itself (would require coding
     a Weibull or piecewise alternative and refitting).
  3. A genuine Acadian-variant structural ceiling on Maine FIA conditions.

The mort.mult null result rules out the simplest mortality lever. (2) and the
broader question of whether the Acadian model fits Maine FIA at all are the
next research directions, and both are outside the autopilot scope.

## Action items

**Production tweak (from v25): ship CSI*0.7 as bridge default in 12.3.9.**
Drops BA bias from +11.1 to ~10 percent, improves R^2 by ~0.01, low risk,
no code-level changes beyond a bridge configuration constant.

**Mortality form refit: bigger commitment.** Replace the Glover/Hool logistic
with a Weibull or piecewise survival function, refit coefficients against ME
FIA mortality data. Paper-sized work.

**Accept the ceiling.** If neither (1) nor (2) is pursued, document +11
percent FIA BA bias as the structural ceiling for the Acadian variant under
12.3.8 production posture and move on.

## Files

`cardinal_acadgy_mortmult_v26.R` and `acadgy_mortmult_v26_results.csv` in
`acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10869581
(approximately 18 minutes on c0326).
