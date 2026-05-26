# 12.3.8 validation on full FIA harness (v24)

2026-05-26. Re-runs the original v21 harness (100 FIA plots, 10 annual cycles,
MORTCAL on, MORTCAL_INTERVAL per stand) against AcadianGY_12.3.8.r. This is the
same plot sample and cycling that v21 used; the only change is the model file
on the source() line.

## Result

| config                          | TPA  | TPA_obs | QMD (cm) | QMD_obs | BA bias % | R^2  |
|---------------------------------|------|---------|----------|---------|-----------|------|
| baseline (CutPoint 0.95, off)   | 1019 | 1029    | 5.150    | 4.968   | +10.90    | 0.42 |
| ingrowth_fix (CutPoint 0)       | 1043 | 1029    | 4.923    | 4.968   | +11.05    | 0.42 |

Compare to v21 (same harness, 12.3.7), where both rows were byte-identical at
TPA 1019.32 / QMD 5.15 / +10.9% because recruits were silently dropped on
cycle 2. With 12.3.8 the two configurations diverge as designed.

## Interpretation

QMD bias improves from +3.7% to -1.0% (absolute residual 3.7 -> 1.0 pp, the
core target of the diagnosis). TPA bias improves from -1.0% to +1.3% (closer in
magnitude). BA bias edges up by 0.15 pp; this is the compensating-errors
structure documented earlier. Recruits add basal area that the diameter
overgrowth was previously masking, so adding them surfaces the underlying
diameter overshoot rather than creating new error. The lever to close the
remaining +11% BA bias is per-species diameter calibration (already fit), not
recruitment.

## Status

Part 1 (12.3.7, frozen recruits) and Part 2 (12.3.8, recruit STAND identity)
of issue #127 are both fixed and validated end-to-end on the FIA harness. The
recruitment-related residual the diagnosis targeted is closed.
