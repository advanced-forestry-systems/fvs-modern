# fvs-modern

fvs-modern is a community maintained fork of the USDA Forest Vegetation Simulator (FVS). It converts the FVS Fortran 77 code base to free-form Fortran 90, builds the regional variants as standalone executables and as shared libraries that Python can drive in process, and carries a Bayesian recalibration of the regional variants against national FIA remeasurement data. A configuration catalog records, per variant, which combination of parameter set and survival model passed long-horizon realism gates.

## Who this site is for

**Users** who want to project stands with an FVS variant from Python, the command line or a web service, choose between default and calibrated parameters, and report the uncertainty of a projection. Start with [Getting started](getting-started/installation.md) and then the [User guide](user-guide/running-projections.md).

**Developers** who build the engine, change the calibration pipeline, regenerate the catalog or cut a release. Start with the [Technical reference](reference/architecture.md).

## Where to start

| If you want to | Read |
|:---|:---|
| Build a library and run one stand in a few minutes | [Quickstart](getting-started/quickstart.md) |
| Decide between default and calibrated parameters, and native or GOMPIT survival | [Choosing a configuration](user-guide/choosing-a-configuration.md) |
| Write keyword files or turn FIA data into stands | [Stand inputs](user-guide/stand-inputs.md) |
| Put a band around a projection | [Uncertainty](user-guide/uncertainty.md) |
| Find out why a run stopped or looks wrong | [Troubleshooting](user-guide/troubleshooting.md) |
| Understand how the pieces fit | [Architecture](reference/architecture.md) |
| Rebuild the configuration catalog | [Validation](reference/validation.md) |
| See how fvs-modern compares with USDA FVS and fvs-conus on the same stands | [Model run-off](reference/runoff.md) |

## The configuration catalog

`fvs2py.catalog.select_config(variant)` reads the catalog in `config/catalog/` and returns the arm to run, the `config_version` to pass to `FVS`, and the environment variables to set. The catalog, its stress test and the per-variant scorecards are archived on Zenodo under the concept DOI [10.5281/zenodo.23777723](https://doi.org/10.5281/zenodo.23777723). The current files are the fifth version (catalog release 1.2.0), regenerated on October 30, 2026; the table is reproduced on [Choosing a configuration](user-guide/choosing-a-configuration.md#current-catalog-by-variant).

## Recent fixes that change results

PR 135, merged on October 30, 2026, fixed four defects in the posterior draw keywords, the in-process calibrated overlay, the FIA stand generator and the draw multiplier basis. Results produced before that date with any of these paths should be regenerated; [Troubleshooting](user-guide/troubleshooting.md#results-made-before-october-30-2026) lists them.

## Reading these pages on GitHub

Every page is plain Markdown under `docs/` and reads on GitHub as well as on the built site. The site navigation is defined in `mkdocs.yml` at the repository root. The [Archive](archive/index.md) holds dated working notes, plans and handoffs that describe how the project got here; they are kept for the record and are not maintained.

## Citing and license

Cite the software through `CITATION.cff` in the repository root. fvs-modern is distributed under the terms in `LICENSE` and `NOTICE`.
