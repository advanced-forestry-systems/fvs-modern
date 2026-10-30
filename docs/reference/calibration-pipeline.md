# Calibration pipeline

The per-variant calibration pipeline fits component models to FIA remeasurement pairs, then writes posterior point estimates into `config/calibrated/<variant>.json` and subsampled posterior draws into `config/calibrated/<variant>_draws.json`. Every script lives in `calibration/R/`, takes `--variant <code>`, and resolves the project root from `FVS_PROJECT_ROOT`, falling back to a path relative to the script.

A second, CONUS-wide pipeline (scripts numbered `30` and up) fits pooled models across all variants and lands its results in the `categories_conus` and `categories_conus_sf` blocks. It is summarized at the end of this page.

## Components fitted

| Step | Script | Component | Method in code | Main outputs in `calibration/output/variants/<variant>/` |
|---|---|---|---|---|
| 01 | `01_fetch_fia_data.R` | Remeasurement pairs | rFIA download, or local FIA tables (`ENTIRE_TREE.csv`, `ENTIRE_PLOT.csv`, `ENTIRE_COND.csv`) | `calibration/data/processed/<variant>/diameter_growth.csv`, `height_growth.csv`, `height_diameter.csv`, `mortality.csv`, `crown_ratio_change.csv` |
| 02 | `02_fit_diameter_growth.R` | Diameter growth, Wykoff ln(DDS) | CmdStan model `calibration/stan/wykoff_dg.stan`; MAP, then HMC, then ADVI | `diameter_growth_samples.rds`, `diameter_growth_posterior.csv`, `diameter_growth_summary.csv`, `diameter_growth_map.csv`, `diameter_growth_method.txt` |
| 03 | `03_fit_height_diameter.R` | Height-diameter, Chapman-Richards `H = 4.5 + a*(1 - exp(-b*DBH))^c` | brms; MAP, then HMC, then meanfield | `height_diameter_samples.rds`, `height_diameter_posterior.csv`, `height_diameter_method.txt` |
| 03b | `03b_fit_height_increment.R` | Height increment (variants with explicit HG parameters: IE, CI, KT, BC, WS, EM) | CmdStan model `calibration/stan/height_increment.stan`; MAP, then HMC, then ADVI | `height_increment_samples.rds`, `height_increment_posterior.csv`, `height_increment_species_index.csv`, `height_increment_method.txt` |
| 04 | `04_fit_mortality.R` | Mortality, logistic | brms; HMC, then meanfield | `mortality_samples.rds`, `mortality_posterior.csv` |
| 05 | `05_fit_crown_ratio.R` | Crown ratio change | brms; HMC, then meanfield | `crown_ratio_samples.rds`, `crown_ratio_posterior.csv` |
| 08 | `08_fetch_stand_data.R` | Stand-level metrics | FIA tree lists | `calibration/data/processed/<variant>/stand_density.csv`, `species_sdi.csv` |
| 09 | `09_fit_stand_density.R` | SDIMAX, BAMAX, self-thinning slope | Quantile regression of ln(TPA) on ln(QMD) at tau 0.95, 0.97 and 0.99, plus a brms hierarchical model | `species_qr_sdimax.csv`, `species_sdimax_calibrated.csv`, `species_bamax_calibrated.csv`, `self_thinning_slopes.csv`, `stand_density_summary.csv`, `stand_density_samples.rds`, `stand_density_overall.csv` |
| 06 | `06_posterior_to_json.R` | Assembly | Reads all of the above | `config/calibrated/<variant>.json`, `config/calibrated/<variant>_draws.json` |
| 07 | `07_diagnostics.R` | Diagnostics | Convergence and prior against posterior | `diagnostics/convergence_summary.csv`, `diagnostics/dg_prior_posterior_values.csv`, plots |

### Inference cascade

The tree-level fits do not all run full HMC. `02_fit_diameter_growth.R` shows the pattern:

1. MAP by L-BFGS (`iter = 10000`), saved to `diameter_growth_map.csv` and used to initialise sampling.
2. HMC with 4 chains, 1,000 warmup and 1,000 sampling iterations, `adapt_delta = 0.95`, `max_treedepth = 12`. The fit is accepted only if fewer than 50 percent of transitions hit maximum treedepth and fewer than 5 percent diverge.
3. Otherwise meanfield ADVI (`iter = 50000`, `output_samples = 4000`).
4. If ADVI also fails, MAP estimates only, with no uncertainty.

The method that produced the fit is written to `<component>_method.txt` (`hmc`, `variational` or `map_only` for diameter growth). Check it before reading convergence statistics: variational fits carry no R-hat or ESS. The brms-based scripts (03, 04, 05) sample 4 chains of 1,500 iterations with 500 warmup before falling back to meanfield; `09_fit_stand_density.R` samples 4 chains of 3,000 with 1,000 warmup.

!!! warning "Unverified"
    `calibration/stan/` (with `wykoff_dg.stan` and `height_increment.stan`) is referenced by the scripts but was not present in the source snapshot used for this page. Model code details beyond what the R scripts show are not documented here.

## Inputs

- FIA tables, either downloaded through rFIA or read from a local directory given by `--fia-dir` or the `FVS_FIA_DATA_DIR` environment variable (scripts 01 and 08).
- The default config `config/<variant>.json`, which supplies the species set, the prior means and the JSON template. `run_variant.sh` exits if it is missing.
- `calibration/data/fvs_variant_ecodivision_map.csv` (read by 01).
- `calibration/data/equation_availability_full.csv`, which gates which components reach `calibration_multipliers`.

## Outputs and how draws are stored

`06_posterior_to_json.R` loads `config/<variant>.json`, writes fitted values into its `categories`, adds a `calibration` metadata block, sources `multipliers.R` to attach `calibration_multipliers`, and writes the result to `config/calibrated/<variant>.json`. If `multipliers.R` fails, the block is skipped with a logged warning and the JSON is still written. See [Configuration files](configuration-files.md) for the schema.

Posterior draws are kept at two levels:

- **Full draws** stay in `calibration/output/variants/<variant>/<component>_samples.rds` (one row per draw, one column per parameter). Script 06 falls back to `<component>_draws.csv` when the RDS is absent.
- **Runtime draws** are written to `config/calibrated/<variant>_draws.json`. Script 06 thins each component to 500 draws by systematic selection (`round(seq(1, n, length.out = 500))`), or keeps all draws when fewer are available, and writes:

```json
{
  "variant": "ne",
  "n_draws": 500,
  "calibration_date": "YYYY-MM-DD",
  "components": ["diameter_growth", "height_diameter", "..."],
  "draws": {"diameter_growth": [{"param": value, "...": "..."}, "..."]}
}
```

Components present can include `diameter_growth`, `height_diameter`, `height_increment`, `mortality`, `crown_ratio` and `stand_density`. `config/uncertainty.py` loads this file (`UncertaintyEngine(variant)`), trims all components to the smallest draw count if they differ, and samples a draw index that it applies to every component, so a single draw keeps parameters from the same row across components.

!!! note "Thinning is per component"
    Each component is thinned independently from its own sample file, and the fits are separate models. Using the same index across components gives a consistent draw label, not a joint posterior across components.

For the CONUS blocks, posterior summaries go into the JSON and full draws stay in CSV files referenced by each block's `draws_csv` field (see `CONUS_BLOCK_SCHEMA.md`).

## Running the pipeline

### One variant

`calibration/slurm/run_variant.sh` runs the full sequence for one variant. It checks that all ten R scripts exist, runs 01 only if `calibration/data/processed/<variant>/` is empty, then runs 02, 03, 03b, 04, 05, 08, 09, 06 and 07 in that order, logging to `calibration/logs/variant_<variant>.log`.

```bash
export FVS_PROJECT_ROOT=/path/to/fvs-modern
bash calibration/slurm/run_variant.sh ne
```

`FVS_PROJECT_ROOT` must be set: the script's fallback is the literal placeholder `/path/to/fvs-modern`. Steps 02 and 06 are fatal on failure, and the script also stops if `diameter_growth_posterior.csv` or the calibrated JSON is missing afterwards. Steps 03, 03b, 04, 05, 08, 09 and 07 log a warning and continue, so a completed run can still lack components; read `components_updated` and `components_skipped` in the `calibration` block.

To rerun a single component and reassemble, call the scripts directly:

```bash
Rscript calibration/R/04_fit_mortality.R --variant ne
Rscript calibration/R/06_posterior_to_json.R --variant ne
Rscript calibration/R/07_diagnostics.R --variant ne
python3 config/validate_calibrated.py
```

Script 06 always starts from the default `config/<variant>.json`, so it rebuilds the calibrated file from every component output present in `calibration/output/variants/<variant>/`, not only the one you refit. It overwrites `config/calibrated/<variant>.json` and `<variant>_draws.json`.

### All variants on SLURM

`calibration/slurm/submit_calibration.sh` is an array job (`--array=1-25`, 4 CPUs, 16 GB, 24 hours) that maps the array index to one of 25 variant codes and calls `run_variant.sh`. Its `--account=default` line must be edited for your allocation, and it loads `R/4.3.0`.

```bash
sbatch calibration/slurm/submit_calibration.sh
```

`calibration/osc/` holds Cardinal-specific submission and resubmission scripts kept for reproducibility of past runs.

## Known gaps

`02_fit_diameter_growth.R` does not write `diameter_growth_species_index.csv`, which `multipliers.R` looks for to map the dense Stan species index onto FVS slots; `03b_fit_height_increment.R` does write its crosswalk. `KNOWN_ISSUES.md` (calibrated configuration defects diagnosed 11 August 2026) traces the unity `dds_multiplier` and `htg_multiplier` arrays in the shipped NE configs to this mapping gap, records that the shipped configs cannot currently be regenerated byte-identically, and notes that the header of `multipliers.R` and `equation_availability_full.csv` disagree on how many variants adopted diameter growth.

## CONUS pipeline

The `30` series builds and fits CONUS-wide models: `30_build_conus_dataset.R` assembles a unified remeasurement dataset with EPA Level 3 ecoregions and site productivity indices; `31_fit_dg_multimodel.R` fits four candidate diameter growth forms; later scripts fit height growth, crown base, mortality (including `34_fit_mortality_gompit.R` and `35_fit_greg_mortality_cch.R`), crown recession and ingrowth. `61_extract_conus_summaries.R` writes summary and draws CSVs under `calibration/output/conus/`, `62_conus_to_variant_json.R` lands per-species `categories_conus.<component>` blocks, and `62b_speciesfree_to_variant_json.R` lands species-free `categories_conus_sf.<component>` blocks. The model forms are specified in `CONUS_MODEL_SPECIFICATION.md` and the JSON blocks in `CONUS_BLOCK_SCHEMA.md`.
