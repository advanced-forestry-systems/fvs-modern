# Uncertainty

The Bayesian calibration produced posterior draws as well as point estimates. fvs-modern propagates them by Monte Carlo: it runs the same stand once per sampled draw, each time with that draw's parameters, and summarises the spread of the resulting trajectories. This page covers what is sampled, how to run an ensemble with `fvs2py`, and what the resulting band does and does not contain.

## Posterior draws

`CALIBRATION.md` describes the draw files: the calibration step `06_posterior_to_json.R` exports 500 posterior draws for each component model to `config/calibrated/<variant>_draws.json`, and all components within one draw share the same MCMC iteration index, which keeps the correlation between, for example, growth and mortality parameters.

`config.uncertainty.UncertaintyEngine(variant, config_dir=None, draws_path=None, seed=None, calibration_dir=None)` loads that file. It looks for `<config_dir>/calibrated/<variant>_draws.json`, or `config/calibrated/<variant>_draws.json` beside the module when `config_dir` is not given; `draws_path` overrides both. The file is read lazily and must carry `n_draws`, `components` and a `draws` object with one list of parameter dictionaries per component. If a component holds fewer draws than `n_draws`, the engine logs a warning and uses the smaller count.

!!! warning "Unverified"
    No `*_draws.json` file is tracked in the repository at commit fb4cf10, so you must produce one with the calibration pipeline (or obtain it from the maintainers) before running an ensemble. The ensemble code below was checked against the source but not executed, and the 500 draw count is taken from `CALIBRATION.md` rather than from a file.

## Running an ensemble

```python
from fvs2py import FVS

fvs = FVS(
    "lib/FVSne.so",
    config_version="calibrated",
    uncertainty=True,
    n_draws=100,
    seed=42,
)
fvs.load_keyfile("runs/quickstart.key")
ensemble = fvs.run_ensemble()      # list of summary DataFrames
band = fvs.uncertainty_summary     # DataFrame, or None if no ensemble

print(band[("atba", "q0500")])                 # posterior median basal area by year
print(band["atba"][["q0025", "q0975"]])        # 95 percent band
```

`run_ensemble(stop_point_code=0, stop_point_year=0)` requires `uncertainty=True`, `config_version` set to `"calibrated"` or `"custom"`, and a keyword file already loaded; it raises `RuntimeError` or `AttributeError` otherwise, and `ImportError` if `config.uncertainty` cannot be imported. For each of `min(n_draws, draws available)` iterations it:

1. samples a draw index uniformly at random, with replacement, from the engine's generator seeded with `seed`;
2. reloads the library and the keyword file;
3. runs to stop point 7 (input read, before imputation);
4. converts the draw to species multipliers on the same basis as the calibrated MAP configuration and sets them with `set_species_attr`: `baimult`, `mortmult` and `htgmult` from `config.calibration_multipliers.multipliers_from_params`, the Python port of `calibration/R/multipliers.R` that also reproduces the shipped MAP multipliers, and `spsdi` as the MAP SDI maximum scaled by the draw's deviation of log SDI maximum from its posterior median; species without a MAP SDI value keep the value FVS already holds;
5. runs the stand to the end with the stop points you passed;
6. stores the summary with two extra columns, `_draw` (iteration number) and `_draw_idx` (posterior draw index).

Because sampling is with replacement, an ensemble of 100 can contain repeated draws. `fvs.uncertainty_results` returns the list from the most recent ensemble.

`compute_multipliers_for_draw` always returns `growth_mult`, `mort_mult` and `hg_mult`, and adds `sdi_values` when the calibrated configuration carries an SDI maximum vector and the draw holds a `stand_density` component. Each multiplier family is gated by `calibration/data/equation_availability_full.csv` and clipped to 0.1 to 10, as in the R serializer. To see them for one draw:

```python
from config.uncertainty import UncertaintyEngine

engine = UncertaintyEngine("ne", seed=42)
draw = engine.get_draw(0)
mults = engine.compute_multipliers_for_draw(draw)   # default_config is no longer used
print(sorted(mults))   # ['growth_mult', 'hg_mult', 'mort_mult', 'sdi_values']
```

!!! success "Fixed in PR 135 (October 30, 2026)"
    Before PR 135, draws were converted by comparing calibrated intercepts with the default FVS coefficients, a different basis from the MAP multipliers, and emitted no `MORTMULT`. The module docstring of `config/calibration_multipliers.py` records that this gave draw `BAIMULT` values of 1.11 to 2.30 on all 108 NE species against MAP values of 0.86 to 1.15 on 58 species. Draws and MAP now share one function and one basis (defect D4). Separately, a draw no longer writes 0 to `spsdi` for species without a value (D2).

## Summarising

`fvs.uncertainty_summary` calls `UncertaintyEngine.summarize_ensemble(results)` with its default quantiles 0.025, 0.25, 0.5, 0.75 and 0.975. The result is indexed by `year` and has two-level columns `(variable, statistic)`. The variables are the [summary columns](outputs.md) in lower case (`tpa`, `atba`, `tcuft` and so on); the statistics are `mean`, `std`, `n_draws` and one column per quantile labelled `q0025`, `q0250`, `q0500`, `q0750`, `q0975`. The `_draw_idx` column is summarised along with the rest; ignore it.

!!! note
    `CALIBRATION.md` shows this table with upper-case variables such as `("BA", "q0500")` and `("VOLUME", ...)`. With `fvs2py` the variables are the lower-case summary column names, so basal area is `("atba", "q0500")`.

For other quantiles, or for a tidy table for plotting:

```python
from config.uncertainty import UncertaintyEngine

band90 = UncertaintyEngine.summarize_ensemble(
    fvs.uncertainty_results, quantiles=[0.05, 0.5, 0.95]
)
long_df = UncertaintyEngine.ensemble_to_long(band90, variables=["atba", "tpa"])
# columns: year, variable, mean, std, n_draws, q0050, q0500, q0950
```

## What the band includes

The band is the spread of the projection across posterior parameter draws, applied as four species-level multipliers on top of the variant's default equations. `CALIBRATION.md` states that it is parametric uncertainty only: it excludes process error (stochastic mortality, regeneration), model structural uncertainty and climate nonstationarity.

It also excludes, by construction:

- variation among stands: the ensemble projects one keyword file;
- uncertainty in the GOMPIT coefficients: `config/greg_mortality_coefficients.csv` holds one row of coefficients per species, so a GOMPIT run carries the same survival model in every draw;
- any component whose multiplier is not produced for your variant (see the diameter growth caveat in [Choosing a configuration](choosing-a-configuration.md)).

The Bakuzis benchmark aggregator (`calibration/python/bakuzis_uncertainty_aggregate.py`) now reports two quantities separately: a within-stand band from each FIA stand's own posterior draws (`bakuzis_uncertainty_stand_bands.csv`, summarised per cell in `bakuzis_uncertainty_cell_bands.csv`) and the between-stand spread of the stands' MAP values. The older band that pooled the draws of all five stands in a cell is still written for backward compatibility, labelled `band_basis = pooled_draws_x_stands`, and the code comments say it mixes parameter uncertainty with stand-to-stand differences. Compare a single-stand `run_ensemble()` band with the within-stand band, not the pooled one.

In the Acadian (ACD) variant the within-stand band at year 100 is near zero for many stands, and that is a property of the engine rather than a sign that the draws were ignored. The mortality multiplier acts only on background mortality, which applies while a stand is below 55 percent of its maximum stand density index (`src-converted/acd/morts.f90` line 522, `IF(RIP .EQ. RN) X=1.0`); above that point density-driven mortality sets the rate. The ACD posterior carries no maximum density term, so every draw approaches the same ceiling, and the draws differ in how fast they get there and in stem counts, not in final basal area. In the October 30, 2026 Paper 1 ensemble 81 of 180 ACD stands gave the same year 100 basal area under all 100 draws while their stem counts still differed in 80 of the 81. In a three stand check a diameter growth multiplier of 1.5 changed year 100 basal area by at most 1 ft² ac⁻¹, while an all-species SDIMAX of 150 or 600 moved all three stands from 178 to 197 ft² ac⁻¹ under the MAP parameters to 84 to 107 or 219 to 278 ft² ac⁻¹ (`validation/bakuzis_stress_v2/d5_20261030/a_acd_band/`).

!!! note "Mortality multipliers and GOMPIT"
    A draw's species `mortmult` does not act on GOMPIT survival. For a species with a GOMPIT coefficient row the mortality routine replaces the background and density rates with the GOMPIT survival term, and the multiplier is not applied (`src-converted/acd/morts.f90` lines 524 to 529, `src-converted/vls/morts.f90` lines 599 to 608). In `vls/morts.f90` a species without a GOMPIT row also loses the multiplier while GOMPIT is on, because it falls back to the background rate without it (line 606). This is read from the source code, not measured in a run.

## Cost

Each draw is a full simulation, so an ensemble of N draws costs about N single runs. `CALIBRATION.md` suggests 50 draws for approximate interactive intervals and 200 to 500 to stabilise the tails, and `calibration/BAKUZIS_UNCERTAINTY_README.md` gives roughly 0.8 to 2.0 seconds per 100 year run on a workstation. The in-process path is faster per draw than the subprocess path used by microfvs.

## The keyword path for draws

`UncertaintyEngine.generate_keywords_for_draw(draw, default_config=None, draw_idx=0, include_comments=True)` writes a draw as `SDIMAX`, `BAMAX`, `MORTMULT`, `BAIMULT` and `HTGMULT` keyword records using the same formatters as the MAP block in `FvsConfigLoader.generate_keywords()`. `BAMAX` is not a posterior quantity and is copied from the MAP configuration, so a draw differs from MAP only in its drawn parameters; `HTGMULT` is written only when some value differs from 1 by more than 0.01. microfvs (`run_fvs(..., uncertainty=True)`) and the Bakuzis benchmark runner use it. In the benchmark runner, posterior draw runs now build FVS with the default parameter set and rely on the draw keyword block; only the MAP arm receives the in-process calibrated overlay, which previously overwrote each draw's SDI maximum with the MAP value.

!!! success "Fixed in PR 135 (October 30, 2026)"
    Before PR 135, `generate_keywords_for_draw` padded every keyword to 16 columns, so FVS read field 1 of `SDIMAX` as blank (all species), each species index landed in the SDI field, and the last record won for every species in every draw (defect D1, per the method's docstring). Draw keyword files written before that date should be regenerated.
