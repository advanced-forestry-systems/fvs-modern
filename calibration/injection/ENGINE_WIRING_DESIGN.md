# Engine wiring design: injecting the species-free equations into FVS

**Date:** 2026-06-14
**Goal:** make the FVS engine run on the trait-driven species-free increment and mortality per tree,
rather than the current per-species multipliers on the native equations. This is the largest remaining
build for Step 3 (the CONUS-wide variant). The fitted forms exist (banked bundles); this document
specifies how they enter the engine and a prototype.

## 1. The mechanism, confirmed in fvs2py

The current calibration path applies FVS keyword multipliers (GROWMULT, MORTMULT, HTGMULT, SDIMAX).
Those act per species per cycle, not per tree, so they cannot carry the within-species, competition- and
trait-driven signal that the species-free equations produce. True per-tree injection needs the FVS API,
which fvs2py already exposes:

- **Stop points.** `FVS.run(stop_point_code=5, stop_point_year=-1)` stops the engine every cycle "after
  growth and mortality have been computed, but prior to applying them" (stop point 5 in the FVS API,
  confirmed in `fvs2py/_base.py`). This is the exact hook: the engine has produced its per-tree
  diameter increment, height increment, and mortality, and they have not yet been applied.
- **Per-tree attribute read/write.** `fvsTreeAttr(name, action, ntrees, values)` (wrapped in
  `fvs2py/_core.py` as `_fvsTreeAttr`) reads or writes per-tree arrays by attribute name (dbh, ht, dg,
  htg, mort/prob, crwdth, species, plot, and the spatial covariates). `action` is "get" or "set".

So the injection loop is: run to stop point 5, GET the tree list and the engine-computed increments,
compute the species-free increments from the banked bundle, SET them back, and resume so FVS applies the
injected values. No Fortran recompile; the native engine still handles bookkeeping, regeneration, and
output.

## 2. The banked bundles (the injection inputs)

Each component is a portable bundle of CSVs plus a manifest (in
`fvs-conus/output/conus/sf_integration/`): `_sf_fixed.csv` (fixed coefficients, posterior means),
`_sf_species.csv` (per-species trait effect on the link scale), `_sf_re_{L1,L2,L3,FT}.csv` (ecoregion
and forest-type random effects), `_sf_gamma.csv` (trait coefficients), and `_sf_manifest.json` (form,
trait columns, covariate names, link). Available: diameter growth (`dg`, `dg_v8_sf`), height growth
(`hg`, `hg_v8rd_sf`), height-to-crown-base (`hcb_v2split`), height-diameter (`htdbh_v2split`), crown
(`cr_t2`), survival (`surv_crz`).

The per-tree prediction is the same linear predictor the offline benchmark uses:

> eta = intercept + trait_effect[species] + z_L1 + z_L2 + z_L3 + z_FT + sum_k(coef_k · covariate_k)

back-transformed by the component link (exp for the log-DDS diameter and height increments, logit for
crown, etc.) and annualized. The covariate construction per component mirrors the `prep_*` functions in
`benchmark_sf_vs_legA.R`, which is the reference implementation.

## 3. The injection loop (per cycle)

1. Initialize FVS for the stand (the existing standinit/treeinit path) with `stop_point_code=5`,
   `stop_point_year=-1`.
2. At each stop:
   a. GET the per-tree arrays needed by the bundles: dbh, ht, cr, species, the competition terms (BAL by
      hardwood/softwood, BA, CCH), and the stand's spatial codes (EPA L1/L2/L3, forest type), plus the
      engine-computed dg and htg.
   b. Compute the species-free dg and htg per tree from the DG and HG bundles, and the survival
      probability from the survival bundle. Annualize to the cycle length.
   c. SET dg and htg to the species-free values; SET the per-tree mortality (prob/tpa) from the survival
      bundle. Optionally blend with the native value using the per-species shrinkage weight
      w = n/(n+kappa) so well-sampled species lean native and rare species lean trait.
   d. Resume the cycle; FVS applies the injected increments and mortality.
3. Continue to the next cycle until the horizon.

This is a predictor-corrector: FVS proposes, the trait model disposes, the engine applies and keeps the
books. The density limit (Step 2) enters here too, as the localized maximum the survival/self-thinning
uses.

## 4. Staged rollout (lowest risk first)

1. **Shadow mode.** Run the loop but do not SET; only log the species-free vs engine increment per tree.
   This validates the predictor against the live engine state with zero behavior change.
2. **Single-component injection.** Inject one component (diameter growth) and benchmark the projected
   stand against observed remeasurement; compare to native and to the multiplier calibration.
3. **Full injection.** Inject diameter growth, height growth, and survival together, with the localized
   maximum, and run the clean benchmark across variants.
4. **Blended injection.** Add the per-species shrinkage blend and tune kappa per component.

Each stage is reversible and benchmarked before the next.

## 5. Why this is the right design

It reaches genuine per-tree, trait-driven behavior without recompiling the FVS Fortran, keeps the native
engine's bookkeeping and output intact, reuses the already-validated bundles and predictor, and is
staged so every step is benchmarked and reversible. It also unifies the two earlier threads: the
localized maximum (Step 2) is the density limit the injected survival uses, and the shrinkage blend is
where the species-dependent and species-free legs meet.

## 6. Validation plan

The offline benchmark already shows the species-free predictions are competitive and well-calibrated
(for example HCB held-out RMSE 0.140 versus 0.147 for the per-species approach, PICP 0.94). The
in-engine validation is the same clean remeasurement benchmark used for the multiplier calibration:
inject, project, compare projected basal area, density, and QMD against observed t2, native, and the
multiplier path. The prototype (`sf_injector.py`) implements the predictor and the stepwise loop; the
one remaining step is an integration run on Cardinal (a compute slot) to execute stages 1 and 2.

## 7. Prototype

`calibration/injection/sf_injector.py` implements: `BundlePredictor` (loads a bundle and computes the
per-tree linear predictor and increment), `SpeciesFreeInjector` (holds the component predictors and
applies them to a tree-list frame), and `run_with_injection` (the stepwise FVS loop using stop point 5
and `fvsTreeAttr`). It is wiring-complete against the confirmed fvs2py API; it needs one integration run
to validate end to end, which is queued for when the cluster clears.
