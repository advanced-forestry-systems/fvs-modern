# Session handoff: AcadianGY 12.3.7 + 12.3.8 ingrowth fix arc

2026-05-27. State at session close.

## What shipped

Two-part fix to a recruitment carry-through bug in the standalone
AcadianGYOneStand path. Both fixes are strict improvements; with INGROWTH = "N"
12.3.8 is byte-identical to 12.3.6.

**12.3.7 (PR #32, merged)** — Recruits in ING.TreeList lacked
dDBH.mult/dHt.mult/mort.mult/max.dbh/max.height columns. Without these,
recruits are frozen at the 3 cm recruitment diameter (dDBH * NA = NA,
coalesces to 0). One block in AcadianGYOneStand before bind_rows: set neutral
multipliers and inherit species size caps from survivors.

**12.3.8 (PR #34, merged)** — Recruits inherit Sum.temp's default
STAND=1/PLOT=1 instead of the parent stand. On any multi-stand harness that
iterates over unique(trees$STAND), the next cycle's per-stand dispatcher finds
no stand_init row for STAND=1, AcadianGYOneStand errors, and the result is
silently NULL'd. Same block: force ingrow$STAND and ingrow$PLOT to the
survivors' single value.

**Direct deployment** — AcadianGY.R in
`ForestVegetationSimulator-Interface-main/fvsOL/inst/extdata/` updated to
12.3.8 with 12.3.7 backed up alongside the earlier 12.3.5/12.3.6 backups.
CHANGELOG entry committed.

## What the validation says

**v24 (FIA, 100 plots, 10 yr, 12.3.8, the headline result)**

| config       | TPA  | QMD cm | BA bias |
|--------------|------|--------|---------|
| baseline     | 1019 | 5.15   | +10.9%  |
| 12.3.8 fix   | 1043 | 4.92   | +11.1%  |
| observed     | 1029 | 4.97   | -       |

QMD bias closes from +3.7% to -1.0%; TPA bias from -1.0% to +1.3%; BA edges
+0.15 pp (recruits surface the diameter overshoot the missing ingrowth was
previously masking).

**v17 (MAGPlot, 262 NB CFI pairs, 10 yr, 12.3.8, 4-cell MORTCAL x CutPoint)**

| config                          | BA bias % | TPH  | QMD cm |
|---------------------------------|-----------|------|--------|
| canonical_off (no MORTCAL, default CP) | -0.04 | 1675 | 12.95 |
| insource_on (MORTCAL on, default CP)   | -6.51 | 1548 | 12.93 |
| insource_on_cp0 (MORTCAL + EV ingrowth)| -6.08 | 1635 | 12.25 |
| **ingrowth_only (no MORTCAL + EV ingrowth)** | **+0.39** | **1764** | **12.30** |

Observed: BA 17.30 m^2/ha, TPH 1807, QMD 11.89 cm.

MORTCAL still over-corrects on Canadian regardless of ingrowth state. The
right operating configuration is MORTCAL off + CutPoint = 0 ingrowth.

**v18 (FIA, 200 plots, 10 yr, 12.3.8, 5-cell lever sweep)**

| config              | BA bias % | TPA  | QMD cm |
|---------------------|-----------|------|--------|
| canonical_off       | +15.35    | 1084 | 5.25   |
| calibrated_off      | +22.26    | 1102 | 5.35   |
| calibrated_on       | +15.16    | 1019 | 5.37   |
| **calibrated_on_cp0**   | **+15.38** | **1042** | **5.16** |
| ingrowth_only       | +15.50    | 1107 | 5.05   |

Observed: BA 98.44 ft^2/ac, TPA 1060, QMD 5.02 in.

Decomposition from canonical_off to production posture calibrated_on_cp0:
calibration +6.91 pp, MORTCAL -7.10 pp, ingrowth +0.22 pp. Within 0.03 pp of
the no-lever baseline. The +15% BA residual is structural and cannot be
closed by these three levers.

## Operating recommendations

**FIA-like Maine**: calibration + MORTCAL + CutPoint = 0 ingrowth.
Best TPA/QMD, BA bias +15% (structural).

**Canadian MAGPlot**: MORTCAL off + CutPoint = 0 ingrowth.
BA bias +0.4%, TPH 1764 (obs 1807), QMD 12.30 cm (obs 11.89).

The split arises because calibration and MORTCAL were fit to Maine
conditions and mechanistically offset there; on Canadian conditions both
levers operate but don't cancel, so the simpler ingrowth-only posture wins.

## What's queued

**Ben Rice** — updated Gmail draft (id r-211259352080030758) with the
comprehensive picture, including both the carry-through fix and the v17/v18
findings. Replace the older draft (id r4590129791148472454) before sending.
Archive `ACD_ingrowth_fix_for_BenRice.zip` (62 KB) on Aaron's desktop;
attach in Gmail before sending.

**Task #131 (stale, in_progress)** — Original "Debug ACD calibrated-path NA
source" from earlier session. Should be reviewable now since the calibrated
path goes through 12.3.8 cleanly in v18.

**Task #180 (pending, new)** — Investigate the structural +15% FIA BA
residual. Three candidates ranked by likelihood:

  1. Climate sensitivity in Kuehne et al. 2020 dDBH equations. The CSI term
     may be under-weighted in the customRun bridge translation. Sensitivity
     analysis: run v25 = v24 with CSI scaled by 0.8 / 1.0 / 1.2 / 1.5 and
     measure BA bias response. If BA closes substantially with higher CSI
     weighting, the bridge is the issue.

  2. Mortality functional form. Current is logistic in size class. Weibull
     or piecewise alternatives may capture the late-stand mortality acceleration
     that the logistic underweights. Would require refitting and is a
     larger commitment than (1).

  3. Accepting +15% as the structural Acadian variant ceiling on Maine FIA.
     The original Kuehne et al. fit was on Acadian-region data including
     New Brunswick CFI; ME FIA is on the edge of that range and may simply
     have systematic differences. If (1) and (2) don't yield substantial
     improvement, document and move on.

## Files of record

`calibration/acd_divergence_2026-05-21/` in fvs-modern main:

- INGROWTH_FIX_12.3.7.md, INGROWTH_FIX_12.3.8.md — diagnosis writeups
- INGROWTH_FIX_12.3.8_v24_validation.md — FIA validation memo
- MAGPLOT_12.3.8_v17_findings.md — Canadian MAGPlot result
- CALIB_12.3.8_v18_findings.md — FIA calibration decomposition
- AcadianGY_12.3.8.r — the model
- patch_ingrowth_fix.py, patch_ingrowth_fix_12_3_8.py — the patchers
- probe_recruit_stand.R — unit test (the strongest evidence)
- cardinal_acadgy_ingrowthfix_v24.R + results CSV — FIA harness
- cardinal_magplot_insource_v17.R + results CSV — MAGPlot harness
- cardinal_acadgy_calib_v18.R + results CSV — FIA calibration harness

PRs: #32 (12.3.7), #34 (12.3.8), #36 (v18 calibration). Direct commits to
main: f0b9817 (v24), 9ad703e (CHANGELOG), 275f610 (v17 MAGPlot). Final
merge commit: 97be982.
