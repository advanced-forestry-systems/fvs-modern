# fvs-modern

[![CI](https://github.com/holoros/fvs-modern/actions/workflows/ci.yml/badge.svg)](https://github.com/holoros/fvs-modern/actions/workflows/ci.yml)
[![Docker](https://github.com/holoros/fvs-modern/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/holoros/fvs-modern/actions/workflows/docker-publish.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19802673.svg)](https://doi.org/10.5281/zenodo.19802673)

A community maintained fork of the USDA Forest Vegetation Simulator, modernized for contemporary computing environments and designed for ease of deployment, variant creation, and calibration.

**Status:** Active development, April 2026
**Origin:** Forked from [USDAForestService/ForestVegetationSimulator](https://github.com/USDAForestService/ForestVegetationSimulator) and [ForestVegetationSimulator-Interface](https://github.com/USDAForestService/ForestVegetationSimulator-Interface)
**Contact:** Aaron Weiskittel, University of Maine (aaron.weiskittel@maine.edu)

## Why this fork exists

FVS is the most widely used individual tree growth and yield model in North America, with a codebase tracing back to 1973. The USDA Forest Service maintains the official version, but the pace of computing modernization has not kept up with the needs of the research and applied forestry community. This fork addresses three specific gaps:

1. **Modern computing standards.** The original codebase is entirely Fortran 77 fixed form with COMMON block global state. This fork converts to free form Fortran 90+, replaces computed GOTOs with SELECT CASE, fixes known runtime bugs, and provides a foundation for eventual module migration.

2. **Turnkey deployment.** Getting FVS-Online running on a new server has historically required deep institutional knowledge. This fork provides Docker containers, AWS AMI builders, platform-specific install scripts (Linux, macOS, Windows/WSL), and automated SSL configuration.

3. **Upstream tracking and variant extensibility.** This fork automatically monitors the USFS GitHub repos for new releases, provides a one-command sync workflow, and includes scaffolding tools for creating new regional variants with calibration templates.

## Repository structure

```
fvs-modern/
    src-converted/              Full FVS codebase in free-form Fortran 90 (2,247 files)
        base/                   Core simulation engine (253 files)
        vbase/                  Virtual base includes (20 files)
        ne/, ie/, cr/, ...      23 US + 2 Canadian regional variants
        fire/, estb/, volume/   Extension modules

    deployment/
        docker/                 Dockerfile, compose, entrypoint (Ubuntu 24.04)
        aws/                    Packer HCL for AMI + EC2 user-data
        scripts/                Platform installers and maintenance tools
            deploy_laptop.sh        Fedora/RHEL deployment
            setup_macos.sh          macOS deployment (Homebrew)
            setup_wsl.sh            Windows/WSL2 deployment
            setup_ssl.sh            Let's Encrypt automation
            check_upstream.sh       USFS GitHub change detection
            sync_upstream.sh        Pull, patch, convert, rebuild, test
            add_variant.sh          New variant scaffolding
            run_regression_tests.sh Full test suite (66/67 passing)
        config/                 NAMESPACE files, DESCRIPTION, YAML config
        microfvs/               REST API integration (FastAPI, based on microFVS)
        fvs2py/                 Python ctypes wrapper (based on Vibrant Planet fvs2py)

    modernization/              Conversion tools and analysis reports
    variant-tools/              Templates for new regional variants
    docs/                       Verification reports, roadmap, update lifecycle
```

## Quick start

For a detailed walkthrough from clone to first projection, see **[docs/getting_started.md](docs/getting_started.md)**.

### Compiler requirement

gfortran is the supported Fortran compiler. Intel ifort is not supported because the `!DEC$ ATTRIBUTES ALIAS` decorators in `src-converted/base/fvs.f90` (and similar locations) emit uppercase symbol exports under ifort that do not match the lowercase + underscore name mangling expected by downstream Fortran callers and Python ctypes wrappers. The build scripts force `FC=gfortran` and will warn if an environment exports anything else.

Install via your platform package manager: `dnf install gcc-gfortran` (Fedora/RHEL), `apt install gfortran` (Debian/Ubuntu), `brew install gcc` (macOS). On HPC systems with environment modules, `module load gcc/12.3.0` works on OSC Cardinal; analogous modules exist on most clusters.

### Option A: Docker (recommended for servers)

```bash
# fvsOL Shiny interface only
cd deployment/docker
docker compose up --build
# FVS-Online available at http://localhost:3838

# Combined: fvsOL + REST API + fvs2py
docker compose -f deployment/microfvs/docker-compose.yml up --build
# FVS-Online at http://localhost:3838
# REST API at http://localhost:8000 (Swagger docs at /docs)
# Jupyter + fvs2py at http://localhost:8888
```

### Option B: Platform-specific install

```bash
# Linux (Fedora/RHEL/Ubuntu)
bash deployment/scripts/deploy_laptop.sh

# macOS
bash deployment/scripts/setup_macos.sh

# Windows (via WSL2)
bash deployment/scripts/setup_wsl.sh
```

### Option C: AWS

```bash
cd deployment/aws
packer build fvs-online.pkr.hcl
# Launch the AMI; FVS-Online auto-starts on port 3838
```

## Upstream synchronization

One of the persistent challenges with FVS forks is keeping up with USFS releases. fvs-modern solves this with automated upstream tracking at three levels:

**Automated weekly checks:** A GitHub Actions workflow (`.github/workflows/upstream-sync.yml`) runs every Monday and compares fvs-modern against the USFS repositories. When new commits are detected, it automatically opens a pull request with the changes for review.

**NVEL volume library tracking:** The upstream [FMSC VolumeLibrary](https://github.com/FMSC-Measurements/VolumeLibrary) is tracked as a git submodule at `upstream/VolumeLibrary/`. An automated audit script numerically compares all 26 coefficient tables against upstream, ignoring F77/F90 formatting differences:

```bash
python3 scripts/audit_nvel_upstream.py
# Current status: 26/26 .inc files match upstream (vollib 20260415)
```

**Manual sync tools:** For immediate updates or when you need to pull specific changes:

```bash
# Check for new upstream commits
./deployment/scripts/check_upstream.sh

# Full sync: pull, patch, convert, rebuild, test
./deployment/scripts/sync_upstream.sh --all
```

See `docs/update_lifecycle.md` for the full update process and `docs/nvel_integration_recommendations.md` for the NVEL tracking strategy.

## Creating a new variant

```bash
./deployment/scripts/add_variant.sh --name me --base ne \
    --description "Acadian variant calibrated for Maine"
```

This scaffolds the directory structure, species mapping template, and calibration parameter files. For runtime calibration without recompilation, use FVS keywords (READCORR, MORTMULT, SITEINDEX) in keyword files.

## Python interfaces

fvs-modern includes two Python interfaces that make FVS accessible without touching Fortran or R.

### fvs2py (ctypes wrapper)

Load any variant as a shared library and drive projections from Python with full access to tree lists, summary tables, and calibrated parameters:

```python
from fvs2py import FVS

fvs = FVS("lib/FVSne.so", config_version="calibrated")
fvs.set_cmdline("--keywordfile=stand.key")
fvs.run()

summary = fvs.summary          # pandas DataFrame: year, TPA, BA, volume
trees   = fvs.tree_attributes  # per-tree DBH, height, crown ratio
```

Uncertainty estimation is built in. Pass `uncertainty=True` to sample from the Bayesian posterior and produce ensemble projections with credible intervals. See `deployment/fvs2py/` for the full API.

### microfvs (REST API)

A FastAPI service that wraps FVS as HTTP endpoints with Swagger documentation:

```bash
cd deployment/microfvs
uvicorn microfvs.main:app --reload
# Swagger docs at http://localhost:8000/docs
```

```bash
# Run a projection via curl
curl -X POST http://localhost:8000/run \
  -H "Content-Type: application/json" \
  -d '{"variant": "ne", "num_cycles": 10, "stand_init": [...]}'
```

The API supports all 25 variants, returns JSON results, and can be deployed behind nginx for multi-user access. See `deployment/microfvs/` for setup.

## Bayesian calibration

The `calibration/` directory contains a complete Bayesian recalibration pipeline for all 25 FVS variants using FIA remeasurement data. Seven model components are refitted per variant: diameter growth (Wykoff), height to diameter allometry (Chapman Richards), height increment (six variants only), mortality (logistic with annualization), crown ratio change, stand density (SDImax), and the self thinning slope (Reineke exponent).

Key features:

- Hierarchical Bayesian models fit via Stan ADVI, producing 500 posterior draws per variant
- Calibrated parameters stored in `config/calibrated/` as JSON (point estimates and full posteriors)
- PERSEUS validation framework comparing projected AGB against observed FIA remeasurements for Maine
- Uncertainty propagation through FVS projections using posterior draw sampling
- Reproducible HPC workflow via SLURM array jobs (see `calibration/slurm/`)

For details on the calibration pipeline, see `CALIBRATION.md`. For the Python projection scripts, see `calibration/python/README.md`.

## Regression test results

| Category | Result |
|----------|--------|
| Shared library builds (25 variants) | 25/25 PASS (23 US + 2 Canadian) |
| API symbol verification | 25/25 PASS (fvssetcmdline_, fvssummary_, fvsdimsizes_, fvstreeattr_) |
| ctypes load under RTLD_LAZY | 25/25 PASS |
| Standalone simulations | 42/42 PASS (iet03 exits cleanly with STOP 10) |
| rFVS API simulation | 1/1 PASS |
| **Total** | **68/68 (100%)** |

## Code modernization status

| Item | Status | Notes |
|------|--------|-------|
| Free form conversion (.f -> .f90) | Complete | 2,247 files |
| FORMAT bug fixes | Complete | 11 instances in vbase/initre.f |
| Computed GOTO -> SELECT CASE | 12/17 in base/vbase | 5 complex dispatchers deferred |
| fvsOL namespace cleanup | Complete | dplyr removed, selective plyr/colourpicker imports |
| COMMON -> module migration | Roadmapped | Phase 1 of modernization plan |

## Related projects

**Official USDA repositories:**

* [USDA FVS](https://github.com/USDAForestService/ForestVegetationSimulator) — original Fortran 77 source, all variants
* [FVS Interface](https://github.com/USDAForestService/ForestVegetationSimulator-Interface) — rFVS and fvsOL (R Shiny) interfaces
* [Open-FVS](https://sourceforge.net/projects/open-fvs/) — earlier community fork on SourceForge (SVN)

**Python interfaces and wrappers:**

* [PyFVS](https://github.com/tharen/PyFVS) — Python package and extension by Tod Haren (USFS), active development
* [PyFVS2](https://github.com/tharen/PyFVS2) — next generation Python interface for FVS (Tod Haren)
* [forest-modeling/PyFVS](https://github.com/forest-modeling/PyFVS) — Python wrappers forked from USFS open-dev branch
* [fvs-python](https://github.com/mihiarc/fvs-python) — pure Python reimplementation of the Southern variant (4 pine species)
* [microFVS](https://github.com/Vibrant-Planet-Open-Science/microfvs) — REST API wrapper (Vibrant Planet)
* [fvs2py](https://github.com/Vibrant-Planet-Open-Science/fvs2py) — Python ctypes wrapper (Vibrant Planet)
* [docker_fvs](https://github.com/Vibrant-Planet-Open-Science/docker_fvs) — Dockerized FVS builds (Vibrant Planet)

**Other forks and tools:**

* [rFVS](https://github.com/SilviaTerra/rFVS) — R interface (SilviaTerra/NCX)
* [open-fvs-mirror](https://github.com/tharen/open-fvs-mirror) — Git mirror of Open-FVS SVN branches
* [UT_FVS](https://github.com/clgiebink/UT_FVS) — tree ring calibration for Utah FVS variant
* [FVSAutomator](https://github.com/mcvittal/FVSAutomator) — batch automation for Ontario variant
* [orenf/fvs](https://github.com/orenf/fvs) — earlier FVS fork (inactive since 2013)

## License

This fork inherits the public domain status of the original USDA Forest Vegetation Simulator. See `LICENSE`.
