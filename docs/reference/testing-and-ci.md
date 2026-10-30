# Testing and CI

This page lists what the GitHub workflows in `.github/workflows/` actually run and the local test commands that exist in the repository. A copy of an older workflow, `modernization/ci/fvs-build-test.yml`, sits outside `.github/` and does not run.

## Workflows

### ci.yml

Triggered on pushes and pull requests to `main`. All jobs run on `ubuntu-latest`.

| Job | What it runs | Gating |
|---|---|---|
| `fortran-syntax` | `gfortran -c -fsyntax-only -ffree-form -fPIC` on six base files: `main.f90`, `myopen.f90`, `volkey.f90`, `algevl.f90`, `evldx.f90`, `errgro.f90` | Fails on any failure |
| `build-ne-variant` | `build_fvs_libraries.sh src-converted lib ne`; checks that `lib/FVSne.so` exists and is at least 1,000,000 bytes; checks `fvssetcmdline_`, `fvssummary_`, `fvsdimsizes_`, `fvstreeattr_` and `volinit_` with `nm`; loads the library through ctypes with `RTLD_LAZY` | Fails if the library is missing or small, or if `fvssetcmdline_` or `fvssummary_` is missing; the other three symbols only warn |
| `build-executables` | Builds NE and LS libraries, then `build_fvs_executables.sh . lib ne ls`, then a banner smoke test that looks for `FVS VARIANT` | `continue-on-error: true`; informational |
| `nvel-audit` | `python3 scripts/audit_nvel_upstream.py --verbose --allow-differ r8clist.inc` with submodules checked out | Fails on unexpected differences |
| `validate-calibration` | `python3 config/validate_calibrated.py` | Fails on any hard error |
| `lint` | `shellcheck --severity=error` on `deployment/scripts/*.sh`; `ruff check calibration/python/ deployment/fvs2py/ config/ --select=E,F,W --ignore=E501` | Neither step fails the job: shellcheck errors are counted and printed, and the ruff command ends in `|| true` |

### cross-platform.yml

Triggered on every push, every pull request and manual dispatch. A matrix over `ubuntu-latest`, `macos-latest` and `windows-latest` builds NE with `build_fvs_libraries.sh src-converted lib ne` and loads the library through ctypes, asserting `fvssetcmdline_`, `fvssummary_` and `fvsdimsizes_`. Linux and macOS install gfortran (Homebrew `gcc@13` on macOS) and load with `RTLD_LAZY`. Windows uses MSYS2 MINGW64, copies the MinGW runtime DLLs next to the library and loads it. Only the Windows leg has `continue-on-error`, although the comment above it says macOS is also non-gating.

The workflow sets `VERBOSE=1` as an environment variable on the build command. The build script ignores that variable and reads only the `--verbose` flag (see [Building the engine](building-the-engine.md#usage)), so these CI logs do not include per-file compile errors.

### docker-publish.yml

Triggered on pushes to `main` that touch `src-converted/**`, `deployment/docker/**`, `deployment/microfvs/**`, `deployment/fvs2py/**`, `deployment/scripts/build_fvs_libraries.sh` or the workflow itself; on published releases; and on manual dispatch. Three jobs build and push images to `ghcr.io` from `deployment/docker/Dockerfile.base`, `deployment/microfvs/Dockerfile` and `deployment/fvs2py/Dockerfile`, named after the repository with no suffix, `-api` and `-fvs2py`. Tags are `latest` on the default branch, `v{{version}}` semver on releases, and the commit SHA.

!!! warning "Unverified"
    The three Dockerfiles named by this workflow were not present in the source snapshot used for this page, so the image contents are not documented here.

### upstream-sync.yml

Runs weekly on Mondays at 12:00 UTC (`0 12 * * 1`) and on manual dispatch. It reads the last-synced dates from `.upstream-sync.json`, counts new commits in `USDAForestService/ForestVegetationSimulator` and `USDAForestService/ForestVegetationSimulator-Interface` since those dates, and when there are any it writes `docs/upstream_report.md`, pushes a branch `upstream-sync/YYYYMMDD` (with an `-HHMM` suffix if that branch exists), and opens a pull request labelled `upstream-sync`. A dependent `build-test` job clones the upstream USFS source, applies the `vbase/initre.f` FORMAT patches, builds `FVSie.so` with upstream's own makefile, installs rFVS, and asserts `maxspecies == 23`. That job tests upstream, not `src-converted/`.

## Local test commands

### Library build and load

```bash
bash deployment/scripts/test_all_variants.sh [SOURCE_DIR] [OUTPUT_DIR] [variant ...]
```

Builds the requested variants (all 25 by default) and loads every library with `RTLD_NOW`, so unresolved symbols fail. Defaults are `src-converted` and `lib`. Its header says the cross-platform matrix uses it, but `cross-platform.yml` does not call it.

### Regression suite

```bash
bash deployment/scripts/run_regression_tests.sh [FVS_BIN_DIR] [TEST_DIR] [VARIANTS...]
```

With no arguments it looks for binaries in `$HOME/fvs-online/bin`, `/srv/fvs/bin` or `$HOME/ForestVegetationSimulator-main/bin`, and for tests in `$HOME/ForestVegetationSimulator-main/tests` or `/srv/fvs/tests`; the per-variant cases in this repository are under `src-converted/tests/`. It runs three parts: rFVS library load tests, standalone simulations of each `*.key` case with the variant executable (120-second timeout per case), and an rFVS API simulation. It exits 1 if any test fails.

`src-converted/tests/test.py` compares two files byte for byte and prints `Comparison Test Passed` or `Comparison Test Failed`.

### Calibrated config validation

```bash
python3 config/validate_calibrated.py
python3 config/validate_calibrated.py --config-dir config --quiet
```

The same check CI runs. See [Configuration files](configuration-files.md#default-and-calibrated) for what it enforces.

### fvs2py unit tests

```bash
cd deployment/fvs2py && pytest
```

The tests in `deployment/fvs2py/fvs2py/tests/` expect a container-style install: `test_base.py` loads `/usr/local/lib/FVSso.so`, and `test_fvs_build.py` runs `/usr/local/bin/FVS<variant>` for every member of `FvsVariant`. `test_core.py` uses mocked libraries. `pyproject.toml` for fvs2py requires Python 3.11 or later.

### Species-free loader test

```bash
python3 config/tests/test_sf_loader.py config/config_loader.py config/calibrated/ne.sf_preview.json
```

Exercises the `conus_sf` and `conus_hybrid` versions of `FvsConfigLoader` against the preview file only.

### Regression tests for the October 30, 2026 fixes

```bash
python3 -m pytest config/tests/test_draws_overlay_standgen.py -q
```

Added in PR 135. Covers the four defects fixed there (D1 draw keyword layout, D2 SDI maximum overlay, D3 condition join in the FIA stand generator, D4 draw multiplier basis) with synthetic fixtures only. It is not run by `ci.yml`.

### NVEL audit

```bash
python3 scripts/audit_nvel_upstream.py --verbose --allow-differ r8clist.inc
```

Compares `src-converted/volume/NVEL` `.inc` files with `upstream/VolumeLibrary` (override with `--local` and `--upstream`). The upstream path is a git submodule; check it out first.

### Syntax check for a changed Fortran file

The pull request template asks for:

```bash
gfortran -fsyntax-only -ffree-form -std=legacy <file>.f90
```

### Lint

```bash
ruff check calibration/python/ deployment/fvs2py/ config/ --select=E,F,W --ignore=E501
shellcheck --severity=error deployment/scripts/*.sh
```

The root `pyproject.toml` configures ruff with line length 100, target `py39`, and rule sets `E`, `F`, `W` and `I`.

### Run guards

`tools/fvs_run_guards.py` validates a keyword file and scans a listing for silent failures such as a `NUMCYCLE` above the compiled `MAXCYC` of 40:

```bash
python3 tools/fvs_run_guards.py --key run.key --out run.out --rc 0
```

Other options are `--numcycle`, `--expected-cycles`, `--summary-rows`, `--db` and `--prgprm`. `tools/fvs_input_guards.py` runs a self-test when executed directly.
