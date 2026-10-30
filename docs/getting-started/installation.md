# Installation

fvs-modern is used from Python through `fvs2py`, a ctypes wrapper that loads one compiled shared library per FVS variant (`FVSne.so`, `FVSacd.so`, and so on). Installing it therefore means three things: building or obtaining the variant libraries, making the `fvs2py` package and the repository's `config/` package importable, and setting the few environment variables that switch optional engine behaviour on.

This page describes the main branch at commit fb4cf10 (October 30, 2026).

## Requirements

| Component | Requirement | Source |
|:---|:---|:---|
| Fortran compiler | `gfortran` only. The build script resets `FC` to `gfortran` and warns if anything else is exported. | `deployment/scripts/build_fvs_libraries.sh` |
| C compiler | `gcc` (taken from `CC`, default `gcc`) | `deployment/scripts/build_fvs_libraries.sh` |
| Python | 3.11 or later, as declared by `requires-python` in the `fvs2py` package | `deployment/fvs2py/pyproject.toml` |
| Python packages | `numpy`, `pandas` at run time; the package's `requirements.txt` also lists `pytest`, `pytest-mock` and `jupyter` | `deployment/fvs2py/requirements.txt` |

The repository README reports Linux as the validated platform, macOS as building and loading from source, and native Windows builds as in progress; Windows users are directed to WSL2 or the Docker image.

!!! note "Python version"
    `docs/getting_started.md` in the repository still says Python 3.9+. The `fvs2py` package metadata requires 3.11, and some modules carry fallbacks for older interpreters (`enums.py` backports `StrEnum`, `common.py` falls back to `typing_extensions` for `ParamSpec`). Use 3.11 or later.

Install the compiler from your platform's package manager: `dnf install gcc-gfortran` (Fedora/RHEL), `apt install gfortran` (Debian/Ubuntu) or `brew install gcc` (macOS). On clusters with environment modules, load a GCC module first; the build script comments name `gcc/12.3.0` on OSC Cardinal.

## Get the source

```bash
git clone https://github.com/advanced-forestry-systems/fvs-modern.git
cd fvs-modern
```

!!! note "Repository address"
    The manuscript draft in the repository gives `advanced-forestry-systems/fvs-modern` as the repository; `CITATION.cff` and the README badges still point at `holoros/fvs-modern`.

## Build the shared libraries

`deployment/scripts/build_fvs_libraries.sh` compiles the free-form Fortran in `src-converted/` into one shared library per variant. Its arguments are the source tree, the output directory and an optional list of lowercase variant codes:

```bash
# one variant
bash deployment/scripts/build_fvs_libraries.sh src-converted lib ne

# several variants
bash deployment/scripts/build_fvs_libraries.sh src-converted lib ne acd ls

# every default variant
bash deployment/scripts/build_fvs_libraries.sh src-converted lib
```

With no variant list the script builds `ak acd bm ca ci cr cs ec em ie kt ls nc ne oc op pn sn so tt ut wc ws bc on`. A `--verbose` flag prints compiler errors. The output directory receives `FVS<variant>.so` for each variant plus two helper libraries, `libfvs_stubs.so` and `libfvs_stubs_final.so`. On Linux each variant library is linked with an `$ORIGIN` run path against those stub libraries, so keep all of them in the same directory.

The script also guards against a build defect that would otherwise produce a library that loads zero trees: it shadows `base/PRGPRM.f90` with the variant's own parameter header during the build and then checks that the `/CONCHR/` common block has the same size in `intree.o` and `blkdat.o`, failing the variant if not. A line such as `[CONCHR] OK: sizeof(COMMON /CONCHR/) matches` in the build log is the expected outcome.

Standalone executables (used by the REST service and other subprocess pipelines) come from a separate script with the same argument pattern. Without a variant list it builds only `ne` and `acd`:

```bash
bash deployment/scripts/build_fvs_executables.sh . lib ne acd
```

See [Building the engine](../reference/building-the-engine.md) for the build internals.

### Check that a library loads

The libraries carry a few unresolved internal NVEL symbols, so they must be opened with lazy binding. `fvs2py` does this itself (`ctypes.CDLL(path, mode=1)`). To check a library by hand:

```bash
python3 -c "
import ctypes, os
so = ctypes.CDLL(os.path.join('lib', 'FVSne.so'), mode=os.RTLD_LAZY)
print(so.fvssetcmdline_, so.fvssummary_, so.fvsdimsizes_)
"
```

On construction `fvs2py` checks for 19 API routines (`fvs`, `fvsAddActivity`, `fvsAddTrees`, `fvsDimSizes`, `fvsEvmonAttr`, `fvsFFEAttrs`, `fvsGetRestartCode`, `fvsGetRtnCode`, `fvsGetICCode`, `fvsSVSDimSizes`, `fvsSetStoppointCodes`, `fvsSetCmdLine`, `fvsSVSObjData`, `fvsSpeciesAttr`, `fvsSpeciesCode`, `fvsStandID`, `fvsSummary`, `fvsTreeAttr`, `fvsUnitConversion`), accepting either the plain name or the lowercase name with a trailing underscore, and raises `ImportError` listing any that are missing.

## Make fvs2py importable

`fvs2py` lives in `deployment/fvs2py/fvs2py/`. The calibrated parameter overlay and the uncertainty engine live in the repository's top-level `config/` directory and are imported as `config.config_loader` and `config.uncertainty`. Both locations must be on the Python path.

The simplest arrangement, and the one used in the [Quickstart](quickstart.md), is to work from the repository root and set `PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/deployment/fvs2py:$PWD"
python3 -c "from fvs2py import FVS; from fvs2py.catalog import select_config; print('ok')"
```

When the repository root is not importable, `FVS.apply_calibrated_config()` catches the `ImportError`, logs `config_loader not found. Install fvs-modern or add config/ to PYTHONPATH.` at warning level and returns an empty dictionary, so the run proceeds with default parameters. `FVS.run_ensemble()` raises `ImportError` in the same situation.

!!! warning "Unverified"
    An editable install (`pip install -e deployment/fvs2py`) was not completed while this page was written: the dependency step, which pulls in `jupyter`, did not finish in the documentation build environment. If you install the package rather than using `PYTHONPATH`, note that `fvs2py.catalog` finds the catalog CSVs by a path relative to its own file (`../../../config/catalog`). That path only resolves inside a repository checkout, so a non-editable install needs `FVS_CATALOG_DIR` (below) or an explicit `catalog_dir=` argument.

## Docker images

The publishing workflow (`.github/workflows/docker-publish.yml`) builds three images to the GitHub Container Registry: a base image with the variant libraries, a REST API image and an `fvs2py` image, named from the repository path with `-api` and `-fvs2py` suffixes for the latter two.

!!! warning "Unverified"
    The Dockerfiles the workflow and `deployment/microfvs/docker-compose.yml` refer to (`deployment/docker/Dockerfile.base`, `deployment/docker/Dockerfile`, `deployment/microfvs/Dockerfile`, `deployment/fvs2py/Dockerfile`) are not in the source snapshot this page was written from, so the image contents, the published image names and the `docker run` commands in the package README could not be checked.

## Environment variables

Every variable below is read somewhere in the code; the file and line are given so you can check behaviour against the source. Variables are grouped by who reads them.

### Read by the engine (Fortran)

The FVS libraries read these with `GETENV` the first time the relevant routine runs in a process. The GOMPIT loader caches its decision (`LDONE` is a `SAVE` variable), so set the variables before the first `run()` in a Python session, not between runs.

| Variable | Effect | Read at |
|:---|:---|:---|
| `FVS_GOMPIT` | Enables the GOMPIT crown-closure survival model in place of native mortality for species with a fitted row. Any non-blank value enables it except one starting with `0`, `n`, `N`, `f` or `F`. | `src-converted/base/gompmort.f90:61` |
| `FVS_GOMPIT_COEF` | Path to the GOMPIT coefficient CSV (columns `SPCD,n,b0,b1,b2,b3,b4`, optional `b5`). Required whenever GOMPIT is on, including when it is switched on with the `GOMPMORT` keyword; if blank or unreadable the engine writes a `GOMPMORT: ... NOT enabled` line to the main output and keeps native mortality. The shipped file is `config/greg_mortality_coefficients.csv`. | `src-converted/base/gompmort.f90:78` |
| `FVS_GOMP_MORTMULT` | Multiplier on the GOMPIT annual hazard (default 1.0; values at or below zero reset to 1.0). | `src-converted/base/gompmort.f90:149` |
| `FVS_GOMP_SFLOOR` | Floor on annual GOMPIT survival for all trees (default 0, off; values outside 0 to 1 reset to 0). | `src-converted/base/gompmort.f90:153` |

GOMPIT is compiled into every variant library except BC (`base/gompmort.f90` is absent from `src-converted/bin/FVSbc_sourceList.txt`).

The engine also reads a set of development switches for the CONUS substitution arms and for the Northeast mortality and height growth diagnostics. All of them default to the frozen behaviour when unset, none is selected by the configuration catalog, and they are listed here only so that a stray value in your shell is recognisable.

| Variable | Compiled into | Read at |
|:---|:---|:---|
| `FVS_GREGDG`, `FVS_GREGDG_COEF` | EM, NE, PN, SN, WC | `src-converted/base/gregdghg.f90:36`, `:52` |
| `FVS_GREG_DD0`, `FVS_GREG_TD`, `FVS_GREG_PPT_SM`, `FVS_GREG_DD18` | EM, NE, PN, SN, WC | `src-converted/base/gregdghg.f90:44` to `:47` |
| `FVS_GREGHG`, `FVS_GREGHG_COEF` | EM, NE, PN, SN, WC | `src-converted/base/greghghg.f90:63`, `:93` |
| `FVS_GREG_EMT`, `FVS_GREG_ELEV`, `FVS_GREG_CONFDIR`, `FVS_GREG_MCW_COEF` | EM, NE, PN, SN, WC | `src-converted/base/greghghg.f90:75`, `:77`, `:80`, `:158` |
| `FVS_GREGCRW`, `FVS_GREGCRW_COEF` | EM, NE, PN, SN, WC | `src-converted/base/greghghg.f90:256`, `:274` |
| `FVS_GREGHTDBH`, `FVS_GREGHTDBH_COEF` | NE | `src-converted/base/greghtdbh.f90:47`, `:55` |
| `FVS_MORT_PMSDIL` | NE | `src-converted/ne/grinit.f90:253` |
| `FVS_HTG_RELHTA_ALPHA`, `FVS_HTG_RELHTA_CAP`, `FVS_HTG_GMODK` | NE | `src-converted/ne/htgf.f90:70`, `:74`, `:78` |
| `FVS_MORT_RELCAP`, `FVS_MORT_VMCRB`, `FVS_MORT_VMCRREF`, `FVS_MORT_VMDIAG` | NE | `src-converted/ne/varmrt.f90:124`, `:130`, `:133`, `:137` |
| `FVS_MORT_HALVE`, `FVS_MORT_CRB`, `FVS_MORT_CRREF` | CS, LS, NE | `src-converted/vls/morts.f90:171`, `:174`, `:177` |

### Read by fvs2py

| Variable | Effect | Read at |
|:---|:---|:---|
| `FVS_CATALOG_DIR` | Directory holding `catalog_recommendation_by_variant.csv` and `catalog_gates_by_arm.csv`, used when no `catalog_dir` argument is passed. Default: `config/catalog` in the checkout. | `deployment/fvs2py/fvs2py/catalog.py:62`, `:75` |
| `FVS_CONFIG_DIR` | Used only by `select_config` to build the `FVS_GOMPIT_COEF` path it returns, when no `config_dir` argument is passed. Default: `config/` in the checkout. | `deployment/fvs2py/fvs2py/catalog.py:109` |

`FVS_CONFIG_DIR` is not read by the `FVS` class or by `FvsConfigLoader`. To point those at a different configuration directory, pass `config_dir=` to the `FVS` constructor.

### Read by the calibration and benchmark drivers only

`FVS_PROJECT_ROOT` and `FVS_LIB_DIR` are read by the driver scripts under `calibration/python/`, not by `fvs2py`. For example, `calibration/python/bakuzis_uncertainty_comparison.py` reads `FVS_PROJECT_ROOT` (default: the repository root, `parents[2]` of the script path) at line 263 and `FVS_LIB_DIR` (default `<root>/lib`) at line 270, and also reads `FVS_CONFIG_DIR` at line 271 and `BAKUZIS_OUTPUT_DIR` at line 273. `calibration/python/fia_stand_generator.py` reads `FIA_DATA_DIR` (default `~/fia_data`) at line 289.

### Read by FVS-Online

`deployment/config/fvsol_config.R` reads `FVSOL_CONFIG_FILE` (line 67) to locate the YAML configuration and lets any `FVSOL_<SECTION>_<KEY>` variable override the corresponding YAML value (line 101).

!!! note "FVS_NOAUTOES is not an environment variable"
    No code reads `FVS_NOAUTOES`. `NOAUTOES` is an FVS keyword (keyword option 72) whose handler (`src-converted/vbase/initre.f90`, line 2815 onward) calls `ESNOAU` to switch off automatic establishment; put it in the keyword file.
