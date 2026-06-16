# fvs-conus: focused status and plan (2026-06-15)

Re-centering on the CONUS unified variant (the species-free and species-dependent component
equations), as distinct from fvs-modern (the per-variant recalibration of the existing engine).
This is "Step 3," the long-term goal: one trait-aware engine for all species, retiring the 20+
regional variants.

## What fvs-conus is

One set of component equations fit on CONUS-wide FIA remeasurement, with two legs blended per species:

- **Species-dependent leg:** per-species parameters where the data support them.
- **Species-independent (species-free) leg:** one trait-and-climate form covering rare and unsampled
  species, where per-species fits run out of data.
- **Blend:** weight w = n/(n+kappa), so well-sampled species keep their own fit and rare species fall
  back to traits.

Components: diameter growth (Kuehne v8), height growth (ORGANON v8rd), height-diameter (Wykoff
v2split), height-to-crown-base, crown ratio (CR2), survival (gompit), ingrowth (count + trait
composition), and the stand-level density layer (built on the localized maximum SDI).

## What is actually done (banked and validated)

- **All six tree-level species-free bundles are fit and banked** in
  `output/conus/sf_integration/`: dg, hg, hcb_v2split, htdbh_v2split, cr_t2, surv_crz (plus alternate
  forms dg_v8_sf, hg_v8rd_sf). Each is a portable bundle (fixed coefficients, per-species trait
  effects, ecoregion and forest-type random effects, trait gammas, manifest).
- **Held-out injection validation (vs the per-species leg-A approach) is now run for three
  components**, on held-out Northeast trees, and the species-free form beats the per-species fallback
  in every one, with well-calibrated prediction intervals, and equals the hybrid (so no per-species
  fallback is even needed):

  | component | species-free RMSE | per-species RMSE | species-free R2 | per-species R2 | PICP-95 |
  |---|---:|---:|---:|---:|---:|
  | Height-to-crown-base | 0.140 | 0.147 | 0.298 | 0.226 | 0.94 |
  | Height growth | 0.315 | 0.318 | 0.114 | 0.094 | 0.93 |
  | Height-diameter | 3.55 | 3.77 | 0.716 | 0.679 | 0.95 |

  This is the core fvs-conus validation: the trait-driven equations generalize as well or better than
  per-species fits on data they have not seen, which is exactly the rare-species case the regional
  variants do not cover. Diameter growth and survival are not yet in the benchmark's predictor
  registry (they need prep functions added) and were validated separately (the Greg head-to-head for
  diameter growth).
- **Diameter-growth head-to-head vs Greg Johnson** (Douglas-fir): species-dependent competitive on
  RMSE; species-free within about 22 percent predicting a species with zero training data.
- **The density layer** uses the localized maximum SDI, validated to predict observed self-thinning
  better than species-weighting.

## What is not done (honest gaps)

- **Engine injection: the blocker is now fixed and the path is open (shadow run is the next step).**
  We found why the injection had never run: fvs2py would not import on the Cardinal Python 3.9
  environment because `_base.py` used Python 3.10+ union type hints (`str | os.PathLike`). A one-line
  lazy-annotations fix (`from __future__ import annotations`) resolves it. Confirmed this session:
  `from fvs2py import FVS` imports, `FVS(lib_path="lib/FVSne.so")` instantiates, `run(stop_point_code,
  stop_point_year)` exposes the stop points, and tree dimensions are accessible. The keyfile is set
  via `fvs.keyfile`. So the mechanism (stop point 5, then per-tree read/write) is now runnable; the
  remaining step is a shadow-mode run (build a one-stand keyfile, stop at point 5, read the tree list,
  log species-free vs engine increment, no SET). The fix is committed on the Cardinal checkout; it
  still needs a clean push to the branch (the branch is ahead of Cardinal's local copy).
- **The unified joint fit is a statistical prototype**, not the production estimation. The current
  bundles are fit per component; the joint estimation of mortality with the localized maximum (so the
  density level is intrinsic) is prototyped, not adopted.
- **Ingrowth species-composition** production fit still needs a compute slot.
- **Uncertainty** (posterior draws to projection-level intervals) is designed; not wired end to end in
  the CONUS engine path.
- **PR #70 (the integration scaffold) is open, not merged.**

## Focused next steps for fvs-conus

1. **Finish the held-out injection validation across all components** (running): one table showing,
   per component, species-free vs hybrid vs per-species on point accuracy and interval calibration.
   This characterizes how good the CONUS equations are and where the trait leg suffices alone.
2. **Shadow-mode injection run** on Cardinal: confirm `sf_injector.py` reads the live tree list and
   the species-free predictions match the bundle, with no behavior change. The first real step toward
   a running unified variant.
3. **Single-component injection benchmark:** inject diameter growth only, project, compare to observed
   and to the native engine; then add height growth and survival.
4. **Adopt the joint mortality-plus-maximum fit** so the density level is intrinsic rather than a
   per-variant scalar.
5. **Finalize ingrowth composition** when a slot frees, and wire uncertainty.

## The honest one-liner

The CONUS species-free equations are fit, banked, and competitive on the one component validated so far
on held-out data; the remaining components are validating now; and the real build that remains is
injecting them into the engine, which is prototyped against the confirmed API but not yet run. That
injection, not more fitting, is the critical path to a running unified CONUS variant.
