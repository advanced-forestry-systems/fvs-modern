# FVS x PERSEUS integration -- handoff

Status: the full-FIADB CONUS FVS projection system is built, validated, and
**live** on the PERSEUS dashboard (https://holoros.github.io/perseus-forest-intelligence/).
This document is the single reference for what exists, how it fits together, the
key findings, and what is left.

**Data deposit (Zenodo):** the processed projection products, pipeline code, and
findings docs are archived and citable.
Concept DOI (latest version): https://doi.org/10.5281/zenodo.20555666 |
v1.0.0: https://doi.org/10.5281/zenodo.20555667 (CC-BY-4.0).
Derived from the fvs-modern software (DOI 10.5281/zenodo.19802672).

## 1. What is live

Three national FVS engines, each with four management scenarios, anchored to FIA
carbon, for all 49 forested states (cls `FVS` on the dashboard):

| engine | model id | character |
|--------|----------|-----------|
| default | `fvs_national_default_v1` | native (Dixon/VARMRT) mortality; over-accumulates with no harvest |
| calibrated | `fvs_national_calibrated_v1` | Bayesian-calibrated growth; 13-32% below default; carries a posterior parameter band on the 7 FIA-anchored states |
| gompit | `fvs_national_gompit_v1` | Johnson national density-dependent mortality; caps and gently declines in late succession |

Scenarios (buckets) per engine: `reserve (no harvest)`, `managed (extensive)`,
`managed (harvest)` (realistic), `managed (intensive)`, plus `harvest_c_yr`.
Metrics: `agc_live_total` (Tg C) and `agb_dry` (Tg). The invariant stress test
(`fvs_dashboard_stress.py`) passes with 0 violations.

The canonical dataset is **v3** (treeinit_h, 100% tree heights). Live commit
chain on `holoros/perseus-forest-intelligence` ends at the v3 refresh.

## 2. The data pipeline (raw FIA -> dashboard)

`run_fvs_perseus_pipeline.sh` chains all of it. Stages:

1. **treeinit fix** (`treeinit_fix_v2.py`): the DataMart `FVS_TREEINIT_PLOT`
   under-expands overstory TPA ~6.5x and is ~30% short on heights. Rejoin each
   row by `TREE_CN` to the raw FIA `<ST>_TREE.csv` to restore `TPA_UNADJ` and fill
   `HT` (FIA modeled height + per-species H-D imputation for the rest). Output
   `FIA_fresh/treeinit_h/`. **Always run this first.** It is the single most
   important fix -- without it the campaign biomass is ~6x too light.
2. **strata inputs**: `build_plantation_flag.py` (STDORGCD -> plt_plantation.csv);
   `build_state_harvest_rates.R` (conus_hcs rasters sampled at plot locations ->
   state_harvest_rates.csv).
3. **campaign** (SLURM arrays, `submit_conus_fvs_v3.slurm`,
   `submit_conus_gompit_v3.slurm`, `--array=0-380%40`): runs every FIA plot
   through the proper regional FVS variant, default + calibrated + gompit arms,
   20x5yr cycles, AGB via NSBE. Output `out_fvs_v3/`, `out_gompit_v3/`
   (`conus_<variant>_b<batch>.csv`). Median batch ~45 min; whole thing overnight.
4. **aggregate** (`fvs_perseus_aggregate.py`): per-state density series
   (Mg/ha), 10/90 plot-percentile band. -> `perseus_series_<cfg>_v3/`.
5. **managed scenarios** (`fvs_managed_v2.py`): plantation-aware. Intensive
   regime only on plantation plots; extensive on natural; calibrated so the area
   blend reproduces the sampled rate. -> `managed3_<cfg>/`.
6. **posterior CI** (`fvs_posterior_uncertainty.py`, SLURM array over
   state x variant x draw; `manifest_post.tsv`): parameter uncertainty per
   anchored state. -> `posterior_ci_all.csv`.
7. **merge + ribbon + stress + push** (on a CLEAN `origin/main` checkout so
   unreleased feature work never leaks): `fvs_perseus_merge.py` (FVS_MANAGED_ROOT
   set) injects engines x scenarios; `fvs_posterior_ribbon.py` applies the
   calibrated parameter band; `fvs_dashboard_stress.py` MUST report 0 before push.

Area model: each state's totals use a fixed area anchored so the 2030 reserve
carbon reproduces `fia.json` tg_agc (7 anchored states: ME GA IN ID MN OR WA;
median ha/plot fallback elsewhere). Series anchored at **2030** because FVS
reports the 2025 inventory treelist before its H-D model runs (the height fix
does not change this -- it is a reporting-order artifact, confirmed).

## 3. Key findings

* **Treeinit expansion bug (the big one):** `FVS_TREEINIT_PLOT` TREE_COUNT
  under-expands overstory ~6.5x, concentrated in eastern variants; caught via a
  QA check before anything went public; fixed against the raw FIA TREE table.
* **gompit caps carbon:** density-dependent national mortality plateaus and
  gently declines late-succession biomass -- the most realistic ceiling of the
  three engines, and why it sits lowest at 2125.
* **Plantation-confined intensive management:** intensive harvest only on FIADB
  plantations (STDORGCD=1, 10.2% of CONUS forest). Where plantations are common
  (GA 28%, OR 20%) the realistic managed path drops below extensive; where
  forests are natural (ME 2.4%) they coincide.
* **Uncertainty hierarchy:** parameter uncertainty (posterior draws) is 0-18% at
  2075, far inside the 30-60% structural spread between engines. The
  mortality-model choice dominates, not the calibrated parameters. The calibrated
  engine is nonetheless informative everywhere (13-32% below default).
* **FIADB vs TreeMap across scales:** the area-expansion choice is negligible at
  CONUS (~2%) but grows with resolution (state 0.36-1.39, forest type 0.23-1.96).
  FVS-on-TreeMap CONUS carbon cross-validates to within 9% of TreeMap's own
  imputed carbon.

## 4. Tooling (calibration/stress/)

| script | role |
|--------|------|
| `treeinit_fix_v2.py` | TPA + height treeinit repair (stage 1) |
| `build_plantation_flag.py` | STDORGCD plantation flag |
| `build_state_harvest_rates.R` | conus_hcs harvest + disturbance rates at plots |
| `run_conus_task_fvstreeinit.py` | per-batch campaign task |
| `fvs_perseus_aggregate.py` | campaign -> per-state density series |
| `fvs_managed_v2.py` | plantation-aware managed scenarios |
| `fvs_perseus_merge.py` | inject engines x scenarios into dashboard series |
| `fvs_posterior_uncertainty.py` | posterior-draw parameter CI (array, --draw-idx) |
| `fvs_posterior_ribbon.py` | apply parameter band to calibrated engine |
| `fvs_strata_trends.py` | landowner/ecoregion/state trends + bootstrap CI |
| `fvs_treemap_conus_compare.py` | FIADB vs TreeMap multi-scale |
| `fvs_dashboard_stress.py` | invariant validator (must be 0 before push) |
| `run_fvs_perseus_pipeline.sh` | orchestrates all stages |

Findings docs: `calibration/output/conus_treeinit_expansion_bug.md`,
`fvs_managed_scenario.md`, `fvs_plantation_scenarios.md`,
`fvs_treemap_conus.md`, `fvs_strata_trends.md`, `fvs_posterior_uncertainty.md`,
`fvs_dashboard_stress.md`, `fvs_perseus_next_steps.md`.

## 5. Cardinal locations

* campaign output: `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_fvs_v3`,
  `out_gompit_v3` (v2 retained for comparison)
* fixed treeinit: `/fs/scratch/.../FIA_fresh/treeinit_h/`
* series / managed / posterior: `perseus_series_*_v3`, `managed3_*`,
  `post_<ST>/`, `posterior_ci_all.csv` under `fvs_stress/`
* dashboard repo: `/fs/scratch/.../perseus_wire` (origin
  `holoros/perseus-forest-intelligence`, deploy from `main` via Vite -> Pages)
* TreeMap VAT areas: `plt_area_treemap.csv`; harvest rates:
  `state_harvest_rates.csv`; plantation flag: `plt_plantation.csv`

## 6. To refresh / reproduce

Run `run_fvs_perseus_pipeline.sh` stage by stage. The push step always uses a
clean `git reset --hard origin/main` so Aaron's `feature/ecoregion-layer` work is
never disturbed; the working tree is restored to that branch afterward. The
stress test gates the push.

## 7. Remaining roadmap (none blocking; all in fvs_perseus_next_steps.md)

* **Dashboard UI cut** to expose the landowner / ecoregion / forest-type
  breakdowns (data + bootstrap CIs already computed; needs a frontend view).
* **HWP carry-over**: route `harvest_c_yr` into a harvested-wood-products pool
  for total-system carbon (the ecoregion-economics layer already models HWP).
* **Stand-age / rotation-aware harvest**: couple the YC owner-rotation logic so
  harvest responds to maturity, and reconcile FVS-managed vs YC-managed.
* **Disturbance temporal basis**: confirm the p_disturbance window (currently
  annualized over 20 yr); fit an annual hazard from the dated layers if possible.
* **Posterior CI for multi-variant western states** (OR/WA span several variants;
  current CI uses the dominant variant only).
* **Climate-sensitive variant** (current runs are climate-static, by agreed scope).
