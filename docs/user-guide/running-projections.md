# Running projections

`fvs2py.FVS` drives one FVS variant library in process through ctypes. A projection is always the same four steps: construct the object on a library, load a keyword file, call `run()`, read results. This page documents the class as it stands at commit fb4cf10 (`deployment/fvs2py/fvs2py/_base.py`). The examples assume the working directory is the repository root and `PYTHONPATH` is set as in [Installation](../getting-started/installation.md).

## The FVS class

```python
FVS(
    lib_path,
    config_version=None,
    config_dir=None,
    custom_config=None,
    uncertainty=False,
    n_draws=100,
    seed=None,
)
```

| Argument | Meaning |
|:---|:---|
| `lib_path` | Path to the variant library, for example `lib/FVSne.so`. |
| `config_version` | Parameter overlay: `None` or `"default"` for the compiled-in parameters, `"calibrated"` for the Bayesian calibrated MAP in `config/calibrated/<variant>.json`, `"custom"` for a JSON file you supply. |
| `config_dir` | Directory holding the variant JSON files and the `calibrated/` subdirectory. Default: the `config/` directory beside `config_loader.py`. |
| `custom_config` | Path to your JSON file; required when `config_version="custom"`. It must follow the schema of `config/<variant>.json`. |
| `uncertainty` | Enables `run_ensemble()` (see [Uncertainty](uncertainty.md)). |
| `n_draws` | Number of posterior draws for `run_ensemble()`. Default 100. |
| `seed` | Seed for the draw sampler. |

The variant code is taken from the library file name: the part of the base name before the first dot, after `FVS`, upper-cased. `lib/FVSne.so` gives `NE`, stored in `fvs.variant`. The configuration overlay uses that code (lower-cased) to find its JSON files, so keep the `FVS<variant>.so` naming.

Construction loads the library and checks that the 19 API routines are present (raising `ImportError` if not). It does not start a simulation.

## Keyword files

`load_keyfile(keywordfile)` reads the file into `fvs.keyfile`, stores its absolute path in `fvs.keyfile_path` and passes `--keywordfile=<absolute path>` to the engine's command line. FVS writes its main output listing beside the keyword file with the extension `.out`.

If the object has already run a simulation, `load_keyfile()` unloads and reloads the library first, so each keyword file starts from a clean engine. It issues a `UserWarning` (`FVS had not completed the previous simulation...`) when the previous run had not reached the finished state.

The engine reads standard FVS keywords; fvs2py does not generate them. [Stand inputs](stand-inputs.md) shows the two ways to supply trees (inline `TREEDATA` records or an SQLite input database), and the [Quickstart](../getting-started/quickstart.md) has a complete inline example.

## Running

```python
fvs.run(stop_point_code=0, stop_point_year=0)
```

`run()` calls the engine repeatedly until it either finishes a stand or reaches a requested stop point. FVS returns control at the end of every stand even when no stop point is set; `restart_code` is then 100. For a keyword file with several stands, call `run()` once per stand. One further call after the last stand lets FVS finish writing the main output; `itrncd` then becomes 2.

When `config_version` is `"calibrated"` or `"custom"`, the first `run()` advances the engine to stop point 7 (input read, before missing values are imputed and before the engine's own calibration), calls `apply_calibrated_config()`, and then continues with the stop points you asked for. You can instead call `apply_calibrated_config()` yourself after stopping at point 7. It returns a dictionary of the species attributes it set (`spsdi`, `baimult`, `mortmult`, `htgmult`), and an empty dictionary when nothing was applied.

!!! warning "Calibrated overlay failures are silent"
    If the configuration file cannot be found or imported, `apply_calibrated_config()` logs a warning through the root logger (for example `Config file not found: .../config/calibrated/ne.json`) and returns `{}`. The run then continues with default parameters and nothing in the summary shows it. Check the return value, or configure logging, before relying on a calibrated run.

!!! warning "The overlay is applied automatically only once per FVS object"
    The flag that records that the overlay was applied is set once and is not reset when `load_keyfile()` reloads the library. A second keyword file loaded into the same `FVS` object therefore runs with default parameters even when `config_version="calibrated"`. In the documentation build, running the quickstart stand twice on one object with a custom configuration gave 172 trees ac⁻¹ and 220 ft² ac⁻¹ in 2125 on the first run and the default 191 and 242 on the second. For every keyword file, stop at point 7 and apply the overlay yourself:

    ```python
    import glob
    from fvs2py import FVS

    fvs = FVS("lib/FVSne.so", config_version="calibrated")
    for key in sorted(glob.glob("runs/*.key")):
        fvs.load_keyfile(key)
        fvs.run(stop_point_code=7)
        applied = fvs.apply_calibrated_config()
        fvs.run()
        print(key, sorted(applied), fvs.summary.iloc[-1]["atba"])
    ```

    With this pattern, and `config_version="custom", custom_config="config/calibrated_sdifix/ne.json"` in place of the calibrated set (which is not in the snapshot used for this page), two successive quickstart runs on one object both gave 172 and 220.

!!! warning "Use one FVS object per library per process"
    ctypes opens a library path once per process and later loads return the same handle, so two `FVS` objects built on the same file share one engine. `load_keyfile()` closes and reopens the library to reset it, but the reset only happens when no other handle holds the library open. In the documentation build a second `FVS` object on `lib/FVSne.so` inherited the first object's finished simulation. Reuse one object for successive keyword files, or run stands in separate processes.

### What the in-process overlay changes

`apply_calibrated_config()` sets at most four species attributes: `spsdi` from the calibrated SDI maximum (species whose calibrated value is missing, `NA` or not positive keep the value FVS already holds), and `baimult`, `mortmult` and `htgmult` computed from differences between calibrated and default coefficient arrays in the JSON `categories` block (`config/config_loader.py`, `apply_to_fvs`). It does not read the `calibration_multipliers` block, which only the keyword path (`FvsConfigLoader.generate_keywords()`) uses. The default NE configuration has no `B0` growth array and no `mortality` category, so for NE the in-process overlay can set only `spsdi`: applying `config/calibrated_sdifix/ne.json` in the documentation build returned `['spsdi']`. Check the returned dictionary for your variant. [Choosing a configuration](choosing-a-configuration.md) compares the two paths.

!!! success "Fixed in PR 135 (October 30, 2026)"
    Before PR 135, `apply_to_fvs` replaced `NA` entries in the calibrated SDI maximum vector with 0 and wrote them, so species without a calibrated value ran with an SDI maximum of 0 (defect D2). `config/catalog/README.md` records that this changed calibrated trajectories in CA, CR, CS, EM, LS and WC, and the catalog's calibrated arms were rerun with the fix (fifth version, PR 137). Calibrated in-process runs made before that date for those variants should be repeated.

### Stop points

Stop points pause the engine inside a cycle so that you can read or change state and then resume with another `run()` call.

| `stop_point_code` | Where FVS stops |
|:---|:---|
| -1 | At every stop location |
| 0 | Never |
| 1 | Just before the first call to the Event Monitor |
| 2 | Just after the first call to the Event Monitor |
| 3 | Just before the second call to the Event Monitor |
| 4 | Just after the second call to the Event Monitor |
| 5 | After growth and mortality are computed, before they are applied |
| 6 | Just before the establishment routines |
| 7 | Just after input is read, before imputation and calibration |

`stop_point_year` is 0 (never), -1 (every cycle, per the `run()` and `set_stop_point_codes()` docstrings; the `stop_point_year` property docstring gives 1 instead) or a calendar year. `set_stop_point_codes(stop_point_code=None, stop_point_year=None)` sets them without running; it raises `ValueError` for a code outside -1 to 7 or a year without a code.

```python
fvs = FVS("lib/FVSne.so")
fvs.load_keyfile("runs/quickstart.key")
fvs.run(stop_point_code=5, stop_point_year=2045)
print(fvs.restart_code)                 # 5
print(fvs.summary["year"].tolist())     # [2025, 2035, 2045]
fvs.run()                               # resume to the end of the stand
```

### Status properties

| Property | Meaning |
|:---|:---|
| `itrncd` | Engine return code: -1 not started, 0 running normally, 1 error (reset by loading new input), 2 all stands finished. |
| `restart_code` | Where the engine last stopped: the stop point code (1 to 7) or 100 after a complete stand. |
| `exit_code` | Worst condition logged so far, from `fvsGetICCode` (0 to 5; see [Troubleshooting](troubleshooting.md)). |
| `stop_point_code`, `stop_point_year` | The values last set, or `None`. |
| `dims` | Dictionary with `ntrees`, `ncycles`, `nplots`, `maxtrees`, `maxspecies`, `maxplots`, `maxcycles`. |
| `stand_ids` | Dictionary with `stand_id`, `stand_cn` and `mgmt_id`; available after a `run()` call. |

## Reading results

`fvs.summary` returns a pandas DataFrame with one row per cycle boundary and the 20 columns described in [Outputs](outputs.md). Rows not yet simulated (all zero) are dropped, and the property returns `None` before any cycle exists. Read it after each stand: the next `run()` or `load_keyfile()` replaces it.

Species information is available through `fvs.species_codes` (columns `fvs_index`, `fvs_alpha`, `fia`, `plants`), `fvs.species_attrs` (one column per species attribute) and `fvs.species`, which joins the two. Individual species attributes are read and written with `get_species_attr(attr)` and `set_species_attr(attr, arr)`, where `arr` must have length `dims["maxspecies"]`. The attribute names are `spccf`, `spsdi`, `spsiteindx`, `bfmind`, `bftopd`, `bfstmp`, `frmcls`, `bfmeth`, `mcmind`, `mctopd`, `mcstmp`, `mcmeth`, `baimult`, `htgmult`, `mortmult`, `mortdia1`, `mortdia2`, `regdmult` and `reghmult`. Set them at stop point 7 or later in the first cycle so the engine uses them:

```python
import numpy as np

fvs = FVS("lib/FVSne.so")
fvs.load_keyfile("runs/quickstart.key")
fvs.run(stop_point_code=7)
mult = np.full(fvs.dims["maxspecies"], 1.2)
fvs.set_species_attr("baimult", mult)
fvs.run()
```

### Tree lists

The class defines `get_tree_attr(attr)`, `set_tree_attr(attr, arr)` and `tree_table(attrs=("id", "species", "dbh", "ht", "dg", "htg", "crwdth"))` for per-tree state, with attribute names `id`, `species`, `dbh`, `dg`, `ht`, `htg`, `crwdth`, `cratio`, `tpa`, `mort`, `plot` and `age`.

!!! warning "Tree attribute access crashes at commit fb4cf10"
    The engine routine `fvsTreeAttr` takes six arguments (`name, nch, action, ntrees, attr, rtnCode`; `src-converted/base/apisubs.f90`, line 129), but `fvs2py` declares and passes five, without `ntrees` (`_base.py`, lines 579 to 597). Every call to `get_tree_attr`, `set_tree_attr` or `tree_table` therefore passes the wrong pointers, and in the documentation build each one ended the Python process with a segmentation fault. Do not call them until the wrapper is fixed.

Until then, write tree lists from the keyword file. The `TREELIST` keyword prints them to the main output listing, and with an SQLite output database the `TREELIDB` keyword inside a `DATABASE` block writes them to a table; the benchmark keyword template in `calibration/python/bakuzis_100yr_comparison.py` uses both.

!!! warning "Unverified"
    Database output (`DATABASE`, `DSNOUT`, `SUMMARY`, `TREELIDB`) was not exercised for this page, because the build used to write it lacked the SQLite C sources.

### Adding trees in memory

`add_trees(dbh, species, ht, cratio, plot, tpa)` passes six equal-length arrays to the engine's `fvsAddTrees` routine and returns the number of records added. Units follow the engine: DBH in inches, height in feet, crown ratio in percent (0 to 100), `species` as the FVS species index, `plot` as a plot number and `tpa` as trees per acre. It raises `RuntimeError` when the engine reports no room or zero trees. Its docstring describes it as a way to load a stand without an input database, called after the engine has read a database-free keyword file. Its argument list matches the Fortran routine, but it was not exercised for this page.

## A complete example

```python
from fvs2py import FVS

fvs = FVS("lib/FVSne.so")
fvs.load_keyfile("runs/quickstart.key")
fvs.run()                         # simulate the stand
summary = fvs.summary.copy()      # keep it before the next run
print(fvs.stand_ids["stand_id"], fvs.exit_code)
fvs.run()                         # let FVS finish its output files
assert fvs.itrncd == 2
```

For many stands, build one keyword file per stand (or one file with several stands and one `run()` per stand), reuse a single `FVS` object, and collect `summary` after each. When you use an overlay, apply it at stop point 7 for each keyword file as shown above.
