# Other interfaces

`fvs2py` is the interface this guide documents in detail. The repository also carries a REST service (microfvs), the R interface used by FVS-Online (rFVS), the FVS-Online web application itself, and an R helper for calibrated keywords. This page says what each one is and where its documentation lives.

## microfvs (REST API)

`deployment/microfvs/` is a FastAPI service based on Vibrant Planet's microFVS. The application object is `microfvs.main:app`; it runs FVS as a subprocess, one stand at a time, using the standalone executables at `/usr/local/bin/FVS<variant>` built by `deployment/scripts/build_fvs_executables.sh`. Inputs go into a temporary SQLite database (`FVS_Data.db`) as `fvs_standinit` and `fvs_treeinit` tables, and the results come back as JSON with the keyword file, return code, standard output, standard error, the output tables and the main listing.

| Method and path | Purpose |
|:---|:---|
| `GET /version` | FVS version running in the service |
| `GET /variants` | Variant codes and names |
| `GET /template` | Default keyword file template (Jinja2) |
| `POST /keyfile` | Render a keyword file from a template and parameters |
| `POST /run` | Run one or more stands from `stand_init` and `tree_init` |
| `POST /outfile` | Return the main output listing for a run |
| `GET /treatments`, `GET /treatments/{treatment_code}` | Treatment keyword components |

The `POST /run` body takes `stand_init`, `tree_init`, `limit`, `template`, `template_params` and `stand_stock_params`. The stand and tree records follow the `FvsStandInit` and `FvsTreeInitRecord` models described in [Stand inputs](stand-inputs.md).

!!! warning "Calibrated parameters are not exposed over HTTP"
    The HTTP `/run` endpoint does not accept `config_version`, `custom_config` or `uncertainty`, and does not pass them on. `docs/getting_started.md` says a `config_version` field is available on the `/run` body; at commit fb4cf10 it is not. The Python function behind the endpoint, `microfvs.utils.run_fvs.run_fvs`, does accept `config_version`, `config_dir`, `custom_config`, `uncertainty`, `n_draws` and `seed`, and with `"calibrated"` or `"custom"` it inserts `FvsConfigLoader.generate_keywords()` before `PROCESS`. Its uncertainty mode uses the draw keyword writer described in [Uncertainty](uncertainty.md#the-keyword-path-for-draws), whose record layout was corrected in PR 135 (October 30, 2026).

The README's example request body (`{"variant": "ne", "num_cycles": 10, "stand_init": [...]}`) does not match the endpoint signature above.

!!! warning "Unverified"
    The keyword file template (`microfvs/templates/keyfile.j2`, loaded by `microfvs.enums.FvsKeyfileTemplate`) and the service Dockerfile are not in the source snapshot used for this page, so the service was not started and the default template's keywords are not documented here.

`deployment/microfvs/docker-compose.yml` defines three services: `fvs-api` (the REST API, host port 8000), `fvs-online` (FVS-Online, port 3838) and `fvs2py` (a Jupyter server with `fvs2py`, port 8888). The Jupyter service is started with `--NotebookApp.token=''`, so it accepts connections without a token; keep that port bound to a trusted network.

## rFVS (R)

rFVS is the R interface from the USDA ForestVegetationSimulator-Interface repository; fvs-modern does not carry its own copy. The deployment scripts install it from a local checkout of that repository (for example `deployment/scripts/deploy_laptop.sh --iface-src /path/to/ForestVegetationSimulator-Interface-main`) and load the fvs-modern shared libraries through it. An R driver that exercises the API ships with the engine tests (`src-converted/tests/APIviaR/Rapi.R`):

```r
library(rFVS)
fvsLoad("FVSie", "../../bin")
fvsSetCmdLine("--keywordfile=base.key")
fvsRun(2, 2030)
fvsGetStandIDs()
```

The same file shows `fvsGetDims()`, `fvsGetTreeAttrs()`, `fvsGetSpeciesCodes()` and `fvsAddActivity()`. Paths in that example are relative to the test directory.

For calibrated parameters from R, `config/fvs_calibrated_keywords.R` defines `load_fvs_config()`, `compute_multipliers()`, `generate_fvs_keywords()`, `inject_keywords_into_keyfile()`, `sample_posterior_draw()`, `generate_draw_keywords()` and `run_uncertainty_ensemble()`.

!!! warning "Unverified"
    The R helper in `config/fvs_calibrated_keywords.R` was not run or reviewed line by line for this page; whether its keyword layout matches the corrected Python writer was not checked.

## FVS-Online

FVS-Online (fvsOL) is the R Shiny web interface. Three guides cover hosting it:

- `deployment/hosting-guide/FVS-Online_Hosting_Guide.md`: production deployment on a public Linux server (Ubuntu 22.04 LTS or RHEL 8/9), with Shiny Server, nginx, TLS and email;
- `deployment/fedora/FVS-Online_Fedora_Deployment_Guide.md`: a local Fedora 39 to 41 deployment, from compiling the libraries through rFVS, fvsOL, the Project Builder and optional nginx, systemd and TLS;
- the README's Docker route (`cd deployment/docker && docker compose up --build`, served on port 3838) and the platform scripts `deployment/scripts/deploy_laptop.sh`, `setup_macos.sh` and `setup_wsl.sh`.

FVS-Online reads its settings from `fvsol_config.yml`, located through `FVSOL_CONFIG_FILE` or a set of default paths, and any value can be overridden with an environment variable named `FVSOL_<SECTION>_<KEY>` (`deployment/config/fvsol_config.R`). The `paths: fvs_bin` setting names the directory of compiled variant libraries.

FVS-Online runs the libraries through rFVS and has no switch for the calibrated overlay or for the configuration catalog; to use either, add the keywords from `generate_calibration_keyfile` to the run's keyword file or set the GOMPIT environment variables for the server process.

!!! warning "Unverified"
    Whether FVS-Online passes environment variables such as `FVS_GOMPIT` through to its R sessions under Shiny Server or the Docker image was not checked.

## Standalone executables

`build_fvs_executables.sh` links each variant with the FVS main program into a standalone executable that reads a keyword file from its command line (`FVSne --keywordfile=stand.key`, the form microfvs and `KNOWN_ISSUES.md` use). Its exit status follows the STOP codes in [Troubleshooting](troubleshooting.md). To apply the calibrated parameters to a standalone run, write the keywords to a file and include them before `PROCESS`:

```python
from config.config_loader import generate_calibration_keyfile

generate_calibration_keyfile("ne", output_path="runs/cal_ne.key")
```

`CALIBRATION.md` suggests including the file with `ADDFILE` or pasting its contents into the keyword file.
