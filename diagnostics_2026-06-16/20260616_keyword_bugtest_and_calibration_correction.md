# Keyword bug test and a correction to the calibration story (2026-06-16)

Aaron was right to be skeptical that widely-used FVS keywords would have "no effect." A proper bug
test at extreme values, with NOCALIB to rule out self-calibration, changes the conclusion.

## Bug test (NE, one stand, 3 cycles, extreme keyword values)

| arm | BA yr0->yr3 | QMD final | verdict |
|---|---|---|---|
| default | 48->140 | 6.73 | baseline |
| NOCALIB | 48->140 | 6.73 | identical -> self-calibration is NOT interfering |
| BAIMULT 0.0 | 48->84 | 5.12 | strong: growth nearly stops. WORKS |
| BAIMULT 2.0 | 48->193 | 8.34 | strong: growth doubles. WORKS |
| CRNMULT 0.2 | 48->140 | 6.73 | ZERO effect |
| CRNMULT 3.0 | 48->140 | 6.73 | ZERO effect -> BUG |
| REGHMULT 0.5 | 48->128 | 6.27 | works |
| REGHMULT 2.0 | 48->160 | 7.38 | works |
| MORTMULT 5.0 | 48->106, TPH 567->400 | -- | works (kills trees) |
| MORTMULT 0.0 | 48->145, TPH 567->589 | -- | works (NE mortality is just low) |

Conclusions:
- BAIMULT, REGHMULT, MORTMULT are FULLY implemented and have strong, proportional effects. My earlier
  "weak leverage" was an artifact of testing mild values (0.75 to 0.95) over a single short interval,
  not a real limitation.
- FVS self-calibration is not damping the multipliers (NOCALIB equals default).
- CRNMULT is genuinely defective in this build: zero effect even at 0.2 and 3.0, though it parses. The
  crown-multiplier arrays exist in CRCOEF.f90 (CRMLT/ICFLG/DHI/DLO) but the wiring from the CRNMULT
  keyword (option 81) into CRCOEF appears incomplete (only blkdat references option 81). This is an
  isolated bug to fix; it matters for crown/volume, less for stand basal area.

## Correction to the earlier calibration conclusion

Earlier I reported that calibration could not improve the engine because the multipliers were inert.
That was half right and half wrong. The EXISTING multipliers are inert because they are relative
(centered on 1) -- that part stands. But the KEYWORDS themselves are powerful, so a properly sized,
SIGNED, per-species/variant calibration does improve predictions. Calibration BAIMULT sweep on NE
(vs observed remeasurement):

| BAIMULT | BA bias | QMD bias |
|---|---|---|
| default | +16.8% | +13.1% |
| 0.5 | +12.0% | +9.9% |
| 0.3 | +10.1% | +8.7% |
| 0.2 | +9.1% | +8.0% |

A signed slowdown around 0.5 to 0.6 removes a third of the basal-area bias and a quarter of the QMD
bias over a single interval, and far more over a full rotation (the bug test: BAIMULT 0.0 takes 30-year
BA from 140 to 84). So Aaron's plan -- calibrate the component equations from FIA by species and
variant and emit as signed multipliers -- is viable through the existing keywords, with no engine
injection required for this (Product A) path.

Note on the residual: even at BAIMULT 0.2 the one-interval BA bias plateaus near +9%, so growth alone
does not zero it on a single interval; the remainder is the standing-stock and likely ingrowth/
accounting terms, and the leverage is much larger over a rotation. Per-species fitting (observed vs
FVS-predicted growth ratio per species) is the refinement on top of the per-variant level.

## fvs2py injection blocker: decisive diagnosis

The FVSne.so contains everything needed -- 248 sqlite symbols (sqlite3_prepare, sqlite3_step, ...),
plus fvs_, filopn_, fvsaddtrees_. So it is NOT a missing-routine or build-stub problem. The standalone
executable (same sources) reads the stand fully; the in-process fvs_ call reads zero trees and writes
empty output (rc=2). The defect is that calling fvs_ via ctypes never runs the DATABASE input phase
that the executable's PROGRAM startup triggers. This needs the fvs2py maintainer's knowledge of the
in-process initialization / stop-point protocol (David Diaz / MicroFVS), or the fvsAddTrees bypass.
It is required only for Product B (full equation replacement), NOT for Product A (the signed-multiplier
calibration), which works now via the verified keywords.

## Revised path

- Product A (deployable now, no injection): calibrate per species and variant from FIA, emit signed
  BAIMULT (and MORTMULT, REGHMULT) -- the keywords work. Fix the CRNMULT wiring so crown calibration is
  also live. This is the immediately viable "calibrated FVS."
- Product B (the fvs-conus equations): still needs the fvs2py init fix to run inside the engine.
