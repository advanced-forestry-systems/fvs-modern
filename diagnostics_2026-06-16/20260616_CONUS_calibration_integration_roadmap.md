# Path to a fully calibrated and integrated CONUS FVS (capstone, 2026-06-16)

The definitive build plan, synthesizing the full verification session. Two products: (A) a calibrated
version of the existing 25 FVS variants via signed component adjustments, and (B) the unified
fvs-conus variant (the new trait-driven equations). Both are gated by specific, identified work; this
document states exactly what is built, what is blocked, and the order.

## Where the truth landed (verified this session)

- FVS over-predicts basal area in 20 of 21 testable variants (+5 to +42%); EC under-predicts. Engine
  benchmark is faithful; the R-equation harness under-predicts because it omits ingrowth/dynamics.
- The existing fvs-modern calibration does not reach or move the engine: multipliers are relative
  (centered on 1, net-zero at stand level); crown and height-diameter were never even emitted; SDIMAX
  was mis-formatted (now fixed/suppressed). So calibrated = default in the engine.
- The refit EQUATIONS are good and validated; they just are not in the engine.
- Uncertainty intervals were parameter-only (cover ~13%); predictive intervals restore ~93% (verified,
  not yet deployed).
- Per-component levers on NE: allometry no stand-BA effect; signed growth (BAIMULT) modest; mortality
  inert on density; the density gap is INGROWTH; localized max SDI is the working density knob.

## Product A: calibrated existing variants via signed, region/species adjustments

Order (dependency-correct, revised by the NE evidence): allometry -> growth -> ingrowth -> max SDI.

| component | keyword / mechanism | status | ceiling |
|---|---|---|---|
| Height-diameter | REGHMULT | format verified; emit signed per species/region | drives volume/biomass, not stand BA |
| Crown ratio | CRNMULT | format verified; emit signed | no stand-BA effect |
| Diameter growth | BAIMULT (signed) | format verified; sized to each variant's over-prediction | WEAK leverage; a scalar cannot fully remove the BA bias |
| Ingrowth / recruitment | ESTAB / FIA ingrowth model | NOT wired into projection; the lever for the density under-prediction | needs wiring (fvs-conus/output/conus/ingrowth fit exists) |
| Mortality | MORTMULT (signed) | format verified | near-inert on NE density; use only where a real signal exists |
| Max SDI | SDIMAX (corrected, localized) | format fixed + verified; localized brms surface exists | co-calibrate level per variant |

Honest expectation: signed adjustments move BA/QMD/TPH modestly in the right direction but cannot fully
calibrate, because the basal-area over-prediction is in the growth EQUATION FORM (a multiplier has weak
leverage) and the density under-prediction is recruitment (ingrowth), which must be added, not
multiplied. Product A is therefore a real but partial improvement; full correction routes through
Product B.

Concrete next builds for A: (1) wire the FIA ingrowth model into the projection and re-test the NE
density gap; (2) fit per-variant signed BAIMULT and localized SDIMAX levels to minimize bias vs
observed; (3) emit CRNMULT/REGHMULT so allometry is live for volume/biomass; (4) re-benchmark all 21
variants. The keyword machinery for all of these is verified and ready.

## Product B: the unified fvs-conus variant (the new equations)

The trait-driven CONUS equations (diameter growth Kuehne v8, height growth ORGANON, height-diameter,
crown, survival, ingrowth) are fit, banked, and validated on held-out FIA. This is the equation set
that actually corrects the bias (it carries the absolute level a multiplier cannot). Integration means
running these equations inside or as the engine.

BLOCKER (maintainer-level): the fvs2py shared-library path does not load a stand. `FVS(lib=FVSne.so)
.load_keyfile().run()` produces zero trees and empty output, while the subprocess executable loads the
same stand fully. Root cause: calling `fvs_` via ctypes runs no input phase at all (empty output,
rc=2), regardless of `--keywordfile` or a stdin-redirect replicating the subprocess; the in-process
entry needs the initialization the executable's program startup performs. Both fix routes exist in the
library (`filopn_`, `fvsaddtrees_`); wiring them is a focused ctypes task that needs the fvs2py
maintainer's knowledge of the in-process init protocol (David Diaz / MicroFVS). Reproductions:
diagnostics_2026-06-16/{clean_repro,fullrun,ctypes_fix}.py.

Until injection works, the fvs-conus equations can only be run in the standalone R/Python projector
(`19_fia_benchmark_engine.R`), whose current harness is stripped (no ingrowth/dynamics) and therefore
under-predicts; it is not yet a faithful CONUS-variant benchmark. So Product B needs either the engine
injection (preferred) or a full-dynamics standalone projector before it can be claimed validated.

## The single ordered plan to "fully calibrated + integrated CONUS FVS"

1. Unblock injection (Product B critical path): fix the fvs2py in-process init so a stand loads
   (focused session with the maintainer). This is the highest-leverage single task; it makes the
   validated equations actually drive predictions.
2. In parallel, finish Product A as the deployable-now improvement: wire ingrowth, fit per-variant
   signed growth + localized max SDI, emit allometry, re-benchmark all 21 variants.
3. Deploy predictive uncertainty intervals (the verified residual-variance fix) so projections carry
   honest coverage.
4. Validate Product B in a full-dynamics framework (post-injection) across all variants vs observed,
   replacing Product A where it wins.

## Status of the machinery (all verified and committed this session)

SDIMAX keyword format fixed; CRNMULT/REGHMULT/BAIMULT/MORTMULT field layouts confirmed; localized
max-SDI surface available; ingrowth model fit; predictive-interval fix verified; engine benchmark
harness and all reproduction scripts committed to `holoros/fvs-modern` branch
`conus-sf-integration-2026-05-21` under `diagnostics_2026-06-16/`. The remaining work is the two builds
above, not more diagnosis.
