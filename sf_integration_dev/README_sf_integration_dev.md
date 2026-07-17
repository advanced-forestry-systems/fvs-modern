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

## conus_greg Python consumer (added 2026-07-17)

`conus_greg_projector.py` is the downstream consumer for the landed Greg CONUS
block (`categories_conus_greg`, read with `FvsConfigLoader(variant,
version="conus_greg")`). It mirrors the conus_sf `project()` entry so both CONUS
options are invoked the same way. Select it and call:

```python
from config.config_loader import FvsConfigLoader
from conus_greg_projector import GregProjector, project
L  = FvsConfigLoader("ne", version="conus_greg")
gp = GregProjector(L)                       # or emt=, td=, emt_td_lookup="stand_emt_td.csv"
gp.dg_annual(97, dbh_in=12, cr=0.5, ht_ft=55, bal=90, elev_ft=1200)   # in/yr, greg_est_dg
gp.hg_annual(97, ht_ft=55, cr=0.5, ccfl=90, cch=0.4, elev_ft=1200)    # ft/yr, greg_est_hg
gp.survival_annual(97, cr=0.5, cch=0.4)                                # annual P, greg_gompit
project(L, "conus_greg", trees0, years=25, step=5)                     # sf-style stepper
```

It evaluates Greg Johnson's deployed forms straight from `species_params`
(DG B0..B6, HG B0=max_ht+B1..B8, survival b0..b4) with the same clamps as the
reference R projector. Species outside Greg's fitted sets (~9-11% of trees) use
the softwood/hardwood median of the in-set coefficients (SPCD<300 = conifer);
pass `fallback_spcd=` to reproduce the reference stand-modal borrow exactly.

CLIMATE (EMT/TD) is required per stand. Supply a scalar, a per-stand CSV
(`STAND_CN,EMT,TD` via `emt_td_lookup=`), or fall back to the documented
ClimateNA 1991-2020 median (EMT=-28.8, TD=24.8 degC). The scalar default is a
placeholder ONLY; production requires a real per-stand lookup (export
`greg_emt_td_lookup.rds` to CSV). `.used_default_climate` flags when the default
was used.

Validation (`test_conus_greg_projector.py`, `demo_conus_greg_ne.py`): on a NE
grid of 5 species x {6,12,20} in DBH, the Python consumer reproduces the R
reference projector (`conus_eq_projector_greg.R`, reading Greg's raw
dg/hg/mort RDS) to max abs diff 1.3e-11 (DG), 7.7e-11 (HG) and 1.5e-6
(survival, bounded by JSON coefficient rounding). DG is positive and monotone
in DBH, species are distinct (Quercus rubra fastest, Acer saccharum slowest).
