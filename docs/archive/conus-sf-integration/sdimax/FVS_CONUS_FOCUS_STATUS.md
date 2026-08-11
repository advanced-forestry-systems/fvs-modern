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

- **Engine injection: environment unblocked, blocker isolated to the fvs2py shared-library path.**
  Update 2 (2026-06-16): the tree-loading issue is now precisely characterized. The exact same stand,
  trees, keyfile, and sqlite database load perfectly through the subprocess executable
  (`run_fvs_projection` auto path: full summary, 639 TPA, 48 ft2/ac BA, SDI 118, projected to 2010 and
  2020). The same inputs through the fvs2py shared-library path (`FVS(lib).load_keyfile().run()`) load
  zero trees (ntrees 0, exit 2). So the data and keyfile are correct; the gap is that the fvs2py
  shared-library path does not initialize FVS to read the stand from the database, even though the .so
  carries the sqlite extension (248 symbols, dbsopen_). That path has never actually been exercised
  (the harness always prefers the subprocess executable), so this is a latent bug in the fvs2py
  wrapper. Two candidate fixes, the next concrete step: (a) supply the FVS input/output unit
  assignments the subprocess provides via stdin (keyfile, tree file, output, treelist) which the
  fvs2py `load_keyfile` does not set, or (b) load the initial tree list via the FVS API `fvsAddTrees`
  instead of the database keyword. Once the .so path loads a stand, stop point 5 plus per-tree
  read/write completes the injection. The stepwise mechanism and the 3.12 environment are proven; this
  one FVS-API initialization gap is the last blocker, and it likely needs a focused session or the
  fvs2py maintainer's input. Diagnostics: `calibration/python/shadow2.py`, `shadow3.py`.

- **(Update 1, 2026-06-15 PM) environment root cause:** the real blocker was the Python environment.
  fvs2py needs Python 3.10+ (it uses union type hints and ParamSpec, 25+ such occurrences); the
  Cardinal default venv is 3.9, so every prior projection silently fell through the harness's
  subprocess fallback, which cannot do stop points or per-tree injection. The `python/3.12` module is
  available; in a clean 3.12 venv (`~/fvs312`, pandas+numpy) fvs2py imports, `FVS(lib_path=FVSne.so)`
  loads the library, and `run(stop_point_code=5, stop_point_year=-1)` executes without error. The one
  remaining issue is that the stand's trees are not reaching the engine through the class path
  (ntrees=0 after the run, exit_code 2 = ran to completion), whereas the same standinit/treeinit load
  fine through the subprocess path. So the next step is a bounded debug: get the class-driven keyfile
  or database load to populate the tree list, then read it at stop point 5. The mechanism and
  environment are proven; this is the last detail before a working shadow run. Shadow test script:
  `calibration/python/shadow_injection_test.py`.

- **(superseded) Engine injection prototype note:** `sf_injector.py` is wiring-complete against the
  confirmed FVS API (stop point 5 plus per-tree attribute read/write).
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
