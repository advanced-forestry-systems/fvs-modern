# Species-free (Leg B) integration: dev/test

Status: DEV ONLY. Nothing here writes the 80 production variant configs.
Hold all production config writes for Aaron's review after the benchmark.

## Engine code
- `config/config_loader.py` (in the repo proper): adds version switch values
  `conus_sf` (species-free, Leg B) and `conus_hybrid` (per species: Leg A where
  a reliable per-species fit exists, else Leg B trait fallback), plus the
  species-free accessors and the runtime trait evaluator
  (`get_conus_sf_runtime_block`, `sf_trait_effect`, `resolve_species_source`,
  `sf_linear_predictor`). Backward compatible; existing versions unchanged.
  Original backed up as `config/config_loader.py.bak_pre_sf`.

## This folder
- `test_sf_loader.py` — unit test for the switch + evaluator. Validated against
  `ne.sf_preview.json`: 16/16 checks pass, evaluator reproduces the precomputed
  trait effect to 1.7e-11.
- `benchmark_sf_vs_legA.R` — held-out three-way benchmark (pure species-free vs
  hybrid vs Leg-A-style) with a consistent base predictor and per-species
  empirical intercepts; reports RMSE/bias/R2 overall, eastern ecoregions, and
  the Leg-A species subset, plus a Leg-A coverage leverage table. Auto-folds new
  bundles via `--scan`. HG predictor registered; CR/HCB/HT-DBH are a small add.
- `61b_extract_speciesfree_summaries.R`, `62b_speciesfree_to_variant_json.R` —
  bundle extractor and the dry-run JSON lander (in calibration/R on Cardinal).

## Hold for review
`62b` defaults to dry-run (writes `{variant}.sf_preview.json`). Do not run it
with `--dry_run=FALSE` against `config/calibrated/` until all component bundles
exist and Aaron approves. Do not merge this branch to main.
