# FVS Parameter Calibration

This document describes the Bayesian calibration pipeline for fvs-modern and
how to use default, calibrated, or custom parameter sets at runtime.

## Overview

FVS ships with parameters estimated from historical data that may not reflect
current growth conditions in a given region. The calibration pipeline fits
Bayesian hierarchical models to FIA remeasurement data, producing updated
parameter estimates with species level partial pooling and full uncertainty
quantification.

Seven component models are calibrated per variant:

| Component | Model | Key Parameters |
|-----------|-------|----------------|
| Diameter growth | Wykoff ln(DDS) | B0, B1 ... B13 per species |
| Height diameter | Chapman Richards | asymptote, rate, shape per species |
| Height increment | Species varying HG | HGLD, HGHC (6 variants only) |
| Mortality | Logistic with annualization | intercept, DBH, BAL, CR effects |
| Crown ratio change | Linear annual delta | BCR1 through BCR4 |
| Stand density (SDIMAX) | Quantile regression + Bayesian | SDICON, R5SDI, or BAMAXA |
| Self thinning slope | Species level Reineke exponent | departure from b = 1.605 |

## Three Parameter Options

### 1. Default (`version="default"`)

Original FVS parameters from Fortran DATA statements. These are the values
that have been in FVS for decades, estimated from data collected primarily
in the 1970s through 1990s.

### 2. Calibrated (`version="calibrated"`)

Bayesian posterior medians fit to national FIA remeasurement data. These
parameters reflect contemporary growth rates, mortality patterns, and
density relationships observed across the most recent FIA measurement
cycles. Produced by the calibration pipeline in `calibration/`.

### 3. Custom (`version="custom"`)

User supplied parameters calibrated to independent data. This option is
designed for users who have access to their own measurement networks and
want FVS to reflect local conditions. Examples include:

  - Cooperative forest inventory networks (CFRU, FIA intensification plots)
  - Long term silvicultural trial data (spacing trials, thinning studies)
  - Regional permanent sample plot networks
  - Company or agency continuous forest inventory (CFI) data
  - Site specific calibrations for management planning

To create a custom calibration, you need a JSON file that follows the same
schema as `config/{variant}.json`. The easiest approach is to copy the
default or calibrated config and modify the parameter arrays you want to
change. See "Creating a Custom Calibration" below.

## Usage

### Python API (fvs2py)

```python
from fvs2py import FVS

# Original FVS parameters (no overlay)
fvs = FVS("path/to/FVSne.so")

# National FIA calibration
fvs = FVS("path/to/FVSne.so", config_version="calibrated")

# Custom calibration from your own data
fvs = FVS("path/to/FVSne.so",
           config_version="custom",
           custom_config="/path/to/my_ne_calibration.json")

fvs.load_keyfile("stand.key")
fvs.run()
```

Parameters are applied automatically at FVS stop point 7 (after input
is read but before imputation and model calibration). This ensures species
attributes like SDI max, growth multipliers, and mortality multipliers are
set before FVS uses them.

### REST API (microfvs)

```python
from microfvs.utils.run_fvs import run_fvs

# Default
result = run_fvs(stand_init, tree_init)

# Calibrated
result = run_fvs(stand_init, tree_init, config_version="calibrated")

# Custom
result = run_fvs(stand_init, tree_init,
                 config_version="custom",
                 custom_config="/path/to/my_ne_calibration.json")
```

When using `calibrated` or `custom`, FVS keywords (SDIMAX, MORTMULT,
GROWMULT, HTGMULT) are injected into the keyfile before the PROCESS
keyword.

### Standalone FVS (Keyfile Approach)

For users running FVS directly from the command line without Python,
generate a keyword snippet and include it in your keyfile:

```python
from config.config_loader import generate_calibration_keyfile

# Write calibration keywords to a file
generate_calibration_keyfile("ne", output_path="cal_ne.key")
```

Then in your main FVS keyfile, add before the PROCESS keyword:

```
ADDFILE
cal_ne.key
```

Or manually append the contents of `cal_ne.key` into your keyfile.

### Comparing Parameter Sets

```python
from config.config_loader import compare_configs

# Print formatted comparison table
print(compare_configs("ne"))
```

This shows default vs. calibrated values for every parameter, with percent
change and number of species affected. Truncated example output for the
Northeast variant:

```
Parameter Comparison: NE (Default vs. Calibrated)
======================================================================
Calibration date: 2026-03-31
Components: diameter_growth, height_diameter, mortality, crown_ratio, sdimax

Parameter                              Default Calibrated     Change
----------------------------------------------------------------------
bark_ratio/BKRAT                          0.91       0.91      +0.0% (0 spp changed)
crown/BCR1                                4.40       4.40      +0.0% (0 spp changed)
growth/B1                                 0.00       0.00      +2.7% (85 spp changed)
growth/B2                                 0.09       0.09      -0.0% (0 spp changed)
other/SIGMAR                              0.05       0.05      +0.0% (0 spp changed)
species_definitions/JTYPE               528.63     528.63      +0.0% (0 spp changed)
...
```

The `(N spp changed)` count reports how many species have a relative
change greater than one percent. For the NE variant, most default
parameters are preserved and the bulk of the calibration effect flows
through the `growth`, `height_diameter`, and `mortality` components
(the last two are applied via SDIMAX, MORTMULT, GROWMULT, and HTGMULT
keyword injection rather than element wise parameter replacement, so
they do not appear as rows in the table).

For programmatic use, `FvsConfigLoader.compare("ne")` returns the same
information as a dictionary, with per parameter arrays of default and
calibrated values, mean percent change, and species change counts.

## Creating a Custom Calibration

### Option A: Modify an Existing Config

1. Copy the default config as a starting point:
   ```
   cp config/ne.json my_ne_calibration.json
   ```

2. Edit the parameter arrays you want to change. For example, to update
   diameter growth intercepts for species 1 through 5:
   ```json
   {
     "categories": {
       "growth": {
         "B0": [-2.1, -1.8, -2.3, -1.9, -2.0, ...]
       }
     }
   }
   ```

3. Use it:
   ```python
   fvs = FVS(lib_path, config_version="custom",
              custom_config="my_ne_calibration.json")
   ```

### Option B: Run the Calibration Pipeline on Your Data

The calibration pipeline in `calibration/R/` can be adapted to fit models
to your own remeasurement data instead of FIA. The key steps are:

1. Format your data to match the schema expected by the fitting scripts
   (see `calibration/data/processed/{variant}/diameter_growth.csv` for
   column definitions)

2. Run the model fitting scripts (02 through 09) pointing at your data

3. Run `06_posterior_to_json.R` to produce the calibrated JSON

4. Use the resulting JSON as a custom config

The required data columns for each component are:

**Diameter growth:** DIA_t1, DIA_t2, years_interval, SPCD, SI, BA, BAL,
CR, SLOPE, ELEV, LAT, LON

**Height diameter:** DIA, HT, SPCD, SI

**Mortality:** DIA, died (0/1), years_interval, SPCD, SI, BA, BAL, CR

**Crown ratio:** CR_init, CR_final, DIA, years_interval, SPCD, SI, BA

**Stand density:** TPA, QMD, BA, leading species, site index

### Option C: Simple Multiplier Approach

If you only need to adjust growth or mortality rates without refitting
the full model, you can create a minimal JSON that contains just the
multiplier keywords. Use `FvsConfigLoader` to generate the keyword block:

```python
from config.config_loader import FvsConfigLoader

loader = FvsConfigLoader("ne", version="calibrated")
keywords = loader.generate_keywords()
# Edit the MORTMULT and GROWMULT values to match your observations
```

## Running the Calibration Pipeline

### On OSC Cardinal (all 25 variants)

Substitute your OSC username and project code below. The example uses
`<OSC_USER>` and `<OSC_PROJECT>` as placeholders; for the maintainer
these are `crsfaaron` and `PUOM0008`, respectively.

```bash
ssh <OSC_USER>@cardinal.osc.edu
cd /users/<OSC_PROJECT>/<OSC_USER>
git clone https://github.com/holoros/fvs-modern.git
bash fvs-modern/calibration/osc/setup_osc.sh     # one time: ~20 min
bash fvs-modern/calibration/osc/submit_cardinal.sh  # submits 25 SLURM jobs
```

Configuration is in `calibration/osc/config_osc.sh`. Adjust
`OSC_ACCOUNT` to your project code and point `FVS_FIA_DATA` (or the
equivalent variable defined in the config) at your FIA database
directory with the expected state subdirectory structure. The same
scripts run on any SLURM scheduler; Cardinal is called out only because
it is the maintainer's environment.

### On a Local Workstation (single variant)

```bash
cd fvs-modern
bash calibration/slurm/run_variant.sh ne
```

### Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 01 | 01_fetch_fia_data.R | Load FIA data, build remeasurement pairs |
| 02 | 02_fit_diameter_growth.R | Bayesian Wykoff model via CmdStan |
| 03 | 03_fit_height_diameter.R | Chapman Richards with species pooling |
| 03b | 03b_fit_height_increment.R | Height growth (6 variants only) |
| 04 | 04_fit_mortality.R | Logistic mortality with annualization |
| 05 | 05_fit_crown_ratio.R | Crown ratio change model |
| 08 | 08_fetch_stand_data.R | Stand level metrics from FIA |
| 09 | 09_fit_stand_density.R | SDIMAX via quantile regression + Bayes |
| 06 | 06_posterior_to_json.R | Merge posteriors into calibrated JSON |
| 07 | 07_diagnostics.R | Convergence checks, prior vs. posterior |

## Output Structure

```
config/
  ne.json                     # Default parameters
  calibrated/
    ne.json                   # Calibrated parameters (from pipeline)
    ne_draws.json             # Full posterior draws (for uncertainty)
  config_loader.py            # Runtime parameter switching module
  uncertainty.py              # Monte Carlo uncertainty propagation engine

calibration/
  output/variants/{variant}/
    diameter_growth_posterior.csv
    height_diameter_posterior.csv
    mortality_posterior.csv
    crown_ratio_posterior.csv
    stand_density_results.rds
    diagnostics/
      convergence_overview.png
      prior_vs_posterior.png
      sdimax_comparison.png
```

## Uncertainty Estimation

A major advantage of the Bayesian calibration is that it produces full
posterior distributions, not just point estimates. This means every
parameter has a distribution reflecting how much the data informed it
relative to the prior. The uncertainty propagation system uses these
distributions to compute credible intervals on FVS projection outputs
(basal area, volume, TPA, mortality, etc.).

### How It Works

The calibration pipeline (06_posterior_to_json.R) exports 500 posterior
draws for each component model into a companion file:
`config/calibrated/{variant}_draws.json`. At runtime, the
`UncertaintyEngine` samples draw indices and assembles complete parameter
vectors, then FVS is run repeatedly with different parameter sets. All
components within a single draw share the same MCMC iteration index, which
preserves correlations between growth and mortality parameters.

The spread of the resulting projection ensemble represents the uncertainty
attributable to parameter estimation. This is parametric uncertainty only;
it does not include process error (stochastic mortality, regeneration),
model structural uncertainty, or climate nonstationarity.

### Usage

#### Python API (fvs2py)

```python
from fvs2py import FVS

# Run 100 posterior draws to get projection confidence bands
fvs = FVS("path/to/FVSne.so",
           config_version="calibrated",
           uncertainty=True,
           n_draws=100,
           seed=42)

fvs.load_keyfile("stand.key")
ensemble = fvs.run_ensemble()

# ensemble is a list of 100 DataFrames (one per draw)

# Get summarized results with credible intervals
summary = fvs.uncertainty_summary
# summary is a pandas DataFrame with a two level MultiIndex on the
# columns: the outer level is the FVS output variable (BA, TPA,
# VOLUME, MORT_BA, etc.) and the inner level is the statistic
# (mean, std, q0025, q0500, q0975). Access a specific band with:
#   summary[("BA", "q0500")]         # posterior median basal area
#   summary[("VOLUME", ["q0025", "q0975"])]  # 95% credible band
```

#### REST API (microfvs)

```python
from microfvs.utils.run_fvs import run_fvs

# Returns a list of FvsResult objects, one per draw
ensemble = run_fvs(stand_init, tree_init,
                   config_version="calibrated",
                   uncertainty=True,
                   n_draws=100,
                   seed=42)
```

#### Manual Draw Selection

For advanced use cases where you want to inspect individual draws or
apply them selectively:

```python
from config.uncertainty import UncertaintyEngine

engine = UncertaintyEngine("ne", seed=42)

# Get a specific draw
draw = engine.get_draw(42)
# draw = {"diameter_growth": {"b0[1]": -2.1, ...}, "mortality": {...}, ...}

# Compute multipliers for this draw
from config.config_loader import FvsConfigLoader
loader = FvsConfigLoader("ne", version="default")
multipliers = engine.compute_multipliers_for_draw(draw, loader.config)
# multipliers = {"growth_mult": array([...]), "mort_mult": array([...]), ...}
```

#### Summarizing Results

```python
from config.uncertainty import UncertaintyEngine

# Summarize an ensemble into credible intervals
summary = UncertaintyEngine.summarize_ensemble(
    ensemble_results,
    quantiles=[0.025, 0.10, 0.50, 0.90, 0.975]
)

# Convert to long format for plotting (e.g., with ggplot2 or matplotlib)
long_df = UncertaintyEngine.ensemble_to_long(
    summary,
    variables=["BA", "TPA", "VOLUME"]
)
```

### Interpreting Uncertainty Bands

Wider credible intervals indicate greater parametric uncertainty. This
typically occurs when:

  1. Fewer FIA observations were available for that species or region
  2. The model structure is a poor fit to the data (high residual variance)
  3. The prior (original FVS parameter) was far from the data, creating
     tension between prior and likelihood

Narrow intervals indicate the data strongly informed the parameter, and
projections are relatively insensitive to remaining uncertainty. For
management planning, the 90% credible interval (5th to 95th percentile)
provides a defensible range of expected outcomes.

### Performance Considerations

Each draw requires a full FVS simulation, so n_draws=100 takes roughly
100x the time of a single run. For interactive use, 50 draws may be
sufficient for approximate intervals. For publication or decision support,
200 to 500 draws are recommended to stabilize the tails.

The fvs2py pathway (shared library) is faster per draw than microfvs
(subprocess) because it avoids process startup overhead. For large
ensembles, fvs2py is strongly preferred.

## Technical Details

The calibration uses Bayesian hierarchical models with species level partial
pooling. This means that species with more FIA observations get estimates
closer to their observed data, while rare species are shrunk toward the
overall mean. This avoids the instability of species specific maximum
likelihood fits while still allowing departure from the pooled estimate
where data support it.

All models are fit using CmdStan (via cmdstanr or brms in R) with 4 MCMC
chains of 2000 iterations each (1000 warmup). Convergence is assessed via
Rhat < 1.01 and bulk ESS > 400 for all parameters.

Priors are informative and centered on the original FVS parameter values,
so the calibrated estimates represent an empirical Bayes update of the
existing FVS parameters rather than fitting from scratch. Where FIA data
are sparse, estimates remain close to defaults.

## Disturbance-aware validation and deployable adjustment levers

See calibration/DISTURBANCE_AWARE_VALIDATION.md for the FIA COND-stratified validation protocol and the four deployable levers (brms site-specific max SDI, sign-aware recruitment injection, signed BAIMULT, and HT-DBH curve correction), with the per-variant config in calibration/validation/data/calibration_config.csv and the all-variant stress test in calib_final.csv. Added 2026-06-17.
