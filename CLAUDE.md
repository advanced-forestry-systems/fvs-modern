# CLAUDE.md: AI Assistant Conventions for fvs-modern

This file orients an AI assistant (or a new contributor) to the repository. For
the user-facing project narrative, the build walkthrough, and the compiler
policy, read `README.md` first; this file covers conventions and the parts of
the tree that the README does not map.

## Project Overview

fvs-modern is a community-maintained fork of the USDA Forest Vegetation
Simulator (FVS), converted from fixed-form Fortran 77 to free-form Fortran 90+
and extended with a Bayesian calibration pipeline that fits regional variants
against national FIA remeasurement data. The project provides:

1. A free-form Fortran 90 conversion of the FVS engine and its regional variants
   (the converted tree lives under `src-converted/`).
2. Deployment infrastructure: Docker, AWS (Packer AMI), and platform install
   scripts for Linux, macOS, and Windows/WSL, plus an automated upstream-sync
   workflow that tracks the USFS repositories.
3. Python and REST wrappers: `fvs2py` (ctypes binding) and `microfvs` (FastAPI).
4. A Bayesian calibration pipeline (Stan + R) producing posterior distributions
   for the component growth/mortality models per variant.
5. An uncertainty-quantification layer that samples posterior draws for ensemble
   projections.

The codebase spans three ecosystems: Fortran (the simulation engine), R (Bayesian
calibration), and Python (APIs, configuration, and post-processing).

## Repository Map

The repository is large (~6,000 tracked files) and has accumulated working
documents at the root. The authoritative source tree is `src-converted/`;
most other top-level items are calibration assets, deployment code, or
project documentation.

```
src-converted/        Converted FVS source. THIS is the build input.
  base/               Core simulation engine
  vbase/              Virtual base includes
  ne/ ie/ cr/ acd/    Regional variants (~23 US + Canadian variants; see README)
    ... and others
  fire/ estb/ volume/ Extension modules
  volume/NVEL/        VENDORED NVEL volume library (mixed .f / .f90 / .c) -- do not "modernize"
  dbsqlite/ dbs/      VENDORED SQLite amalgamation (sqlite3.c/.h) -- do not edit by hand
  archive/            Retired variants and experiments (still fixed-form .f in places)
  stubs/              Generated stub sources used by the build
  tests/              Regression tests (test.py, comparison scripts)

deployment/
  scripts/            Build + install automation
    build_fvs_libraries.sh   Compile variants to .so shared libraries
    run_regression_tests.sh  Full regression suite
    add_variant.sh           Scaffold a new variant
    setup_macos.sh / setup_wsl.sh / deploy_laptop.sh   Platform installers
    setup_ssl.sh, check_upstream.sh, sync_upstream.sh
  docker/             Dockerfile + compose (Ubuntu 24.04, Shiny + R)
  aws/                Packer HCL for AMI + EC2 user-data
  fvs2py/             Python ctypes wrapper
  microfvs/           FastAPI REST interface
  patches/            Auto-applied diffs for upstream compatibility
  config/, services/, fedora/, hosting-guide/   Install support

calibration/
  R/                  Bayesian fitting scripts
  stan/               Stan model specifications (diameter growth, mortality, ...)
  python/             Projection engines and aggregators
  osc/, slurm/        SLURM submission templates for OSC Cardinal (historical)
  data/               FIA remeasurement pairs, processed by component
  output/             Posterior draws, diagnostics, figures, benchmark reports
  analysis/, stress/, figshare/   Analysis notebooks, stress tests, archive bundles

config/
  *.json              Default parameters per variant
  calibrated/         Posterior medians for calibrated variants
  config_loader.py    Runtime parameter switching (default/calibrated/custom)
  uncertainty.py      Posterior-draw sampling and aggregation
  fvs-modern.env.example   Template for the FVS_* environment variables

docs/                 Getting-started, verification reports, roadmap, lifecycle
variant-tools/        Templates and helpers for new variants
modernization/        Conversion tooling and conversion reports
```

### Root-level documents

Several Markdown design/handoff documents currently live at the repository root
rather than under `docs/`: `CALIBRATION.md`, `KNOWN_ISSUES.md`, `STOP_CODES.md`,
`CHANGELOG.md`, and the CONUS-variant set (`CONUS_BLOCK_SCHEMA.md`,
`CONUS_DEVELOPMENT_NOTES.md`, `CONUS_MODEL_SPECIFICATION.md`,
`FVS_CONUS_INTEGRATION_PLAN.md`, `HANDOFF_CONUS_PROJECTION.md`). Treat these as
authoritative for their topics. Consolidating them into `docs/` is a pending
cleanup; until then, check the root before assuming a topic is undocumented.

Note: the working tree may also contain untracked scratch (build artifacts such
as `*.mod`, generated `perseus_*` figures, a local `package.json` / `node_modules`
for an ad-hoc document generator). None of that is tracked or part of the build --
ignore it, and never `git add -A` blindly.

## Build Instructions

Build all variant shared libraries from the converted tree:

```bash
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib
```

The script takes SOURCE_DIR and OUTPUT_DIR as required arguments and produces
per-variant `.so` files (`FVSne.so`, `FVSie.so`, ...) in `./lib`. It checks
dependencies, compiles `base/` to objects, links variant modules against base,
generates the required stub libraries, and builds position-independent shared
libraries for ctypes/fvs2py.

Subset builds:

```bash
bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib ne ak ie
```

**Compiler policy: gfortran only.** Intel `ifort` is not supported. The
`!DEC$ ATTRIBUTES ALIAS` decorators (e.g. in `src-converted/base/fvs.f90`) emit
uppercase symbol exports under ifort that do not match the lowercase-plus-
underscore name mangling expected by downstream Fortran callers and the Python
ctypes wrappers. The build scripts force `FC=gfortran` and warn on anything else.

## Key Conventions

**Fortran:** The converted FVS engine under `src-converted/` is free-form
Fortran 90 (column-independent, `!` comments, `SELECT CASE` in place of computed
GOTOs, modern declarations like `integer, parameter :: MAX_TREES = 10000`).
Exceptions you must not "modernize": the vendored SQLite amalgamation
(`dbsqlite/sqlite3.c`), the vendored NVEL volume library (`volume/NVEL/`, which
retains fixed-form `.f` and C sources), and retired sources under
`src-converted/archive/`. Some `rd/` files are INCLUDE-only declaration units
despite a `.f90` extension -- the build script deliberately skips compiling them
as standalone units (see `KNOWN_ISSUES.md`).

**Python 3.9+:** Type hints preferred for new code. Use `pathlib.Path` for path
operations and `os.environ.get()` with sensible defaults for configuration.

**R (tidyverse):** Pipe operators, tidy data frames (one observation per row),
`readr` for I/O. Avoid global assignment except for data loading.

**No hardcoded HPC paths in new code.** Cardinal paths (`/users/PUOM0008/`,
`/scratch/`) are confined to `calibration/osc/` and `calibration/slurm/` SLURM
templates, kept for reproducibility of past runs. Production code reads
`FVS_PROJECT_ROOT`, `FVS_LIB_DIR`, and `FVS_FIA_DATA_DIR` from the environment
(see `config/fvs-modern.env.example`). Do not extend the hardcoded-path pattern.

**Parameter versioning.** FVS instances accept `version="default|calibrated|custom"`.
Runtime switching goes through `config/config_loader.py` so parameter sources stay
consistent. The uncertainty layer (`config/uncertainty.py`) samples complete
posterior parameter vectors per draw, preserving within-draw parameter
correlations while allowing ensemble spread across draws.

## Testing

```bash
bash deployment/scripts/run_regression_tests.sh
```

This runs library-load tests (ctypes import in Python, `gctorture` in R),
standalone keyword-file simulations, rFVS `.Fortran()` API tests, and
comparative benchmarks (default vs. calibrated parameters). For the live
pass/fail count, rely on the CI badge in `README.md` rather than a number
hardcoded here -- the suite changes as variants are added. The one long-standing
failure is the `FVSie` `iet03` standalone segfault (tracked as issue #5).

Python unit tests:

```bash
cd deployment/fvs2py && pytest
```

## Working With This Repo

- The default branch is `main`. Feature work happens on topic branches; an
  automated workflow opens dated `upstream-sync/YYYYMMDD` PRs when the USFS
  repositories change -- those are bot PRs, not human work.
- Do not commit large binaries (`.docx`, `.pptx`, `.xlsx`, `.pdf`, large `.png`,
  data dumps). Manuscripts, slide decks, and generated figures belong in a
  release asset, Zenodo deposit, or Git LFS -- not in version history. (See the
  history-size note in `docs/HISTORY_REWRITE_PLAN.md`.)
- Releases are tagged `vYYYY.MM.N` with notes; `CHANGELOG.md` is maintained.

## Licensing

The original FVS Fortran source is a USDA Forest Service work and is in the U.S.
public domain (17 U.S.C. section 105); see `NOTICE`. New contributions in this
repository -- the calibration pipeline, analysis and visualization scripts, CONUS
variant equations, and associated documentation -- are MIT licensed; see
`LICENSE`.
