# Stand inputs

FVS reads a stand from its keyword file. The trees can come either from fixed-format records inside the keyword file (`TREEDATA`) or from an SQLite input database with the tables `FVS_StandInit` and `FVS_TreeInit`. Both routes work with `fvs2py`; the repository uses the database route for its Northeast and Acadian benchmark runs and the inline route for other variants (`calibration/python/bakuzis_uncertainty_comparison.py`, `run_scenario`). All tree and stand values are in the engine's imperial units: DBH in inches, heights in feet, expansion in trees per acre.

## Route 1: inline tree records

`calibration/python/inventory_keyfile.py` builds a complete keyword file from two pandas DataFrames. It is the simplest way to produce valid fixed-column records.

```python
import sys
sys.path.insert(0, "calibration/python")

import pandas as pd
from inventory_keyfile import make_inventory_keyfile

stand = pd.DataFrame([{
    "site_index": 60, "age": 45, "aspect": 0, "slope": 10,
    "elevft": 800, "forest_type": 0,
}])
trees = pd.DataFrame({
    "tree_count": [120, 80, 60, 40, 30, 20],   # trees per acre per record
    "species":    [12, 97, 316, 371, 129, 12],  # FIA species codes
    "diameter":   [5.2, 7.8, 6.5, 9.1, 12.4, 3.1],  # in
    "ht":         [38, 52, 45, 58, 70, 24],     # ft
    "crratio":    [55, 45, 40, 45, 50, 60],     # percent
})
key = make_inventory_keyfile(stand, trees, "QS001", "ne", inv_year=2025, num_cycles=20)
open("runs/gen.key", "w").write(key)
```

`make_inventory_keyfile(stand_df, tree_df, stand_id, variant, inv_year=2000, num_cycles=20, calibration_keywords="")` writes:

- `DESIGN` with a basal area factor of -1, one 1 acre fixed plot and a break DBH of 999, so each record's expansion is its `tree_count` in trees per acre;
- `STDINFO` with a per-variant forest code (`STDINFO_DEFAULTS`: 919 for NE and ACD, 904 for LS, 612 for PN, 803 for SN, 110 for IE; other variants fall back to the PN entry), stand age, aspect, slope and elevation in hundreds of feet;
- `SITECODE` with the first tree's species as site species and `site_index` as its site index;
- `INVYEAR`, `NUMCYCLE` and a fixed `TIMEINT` of 5 years, so `num_cycles=20` is a 100 year run;
- any `calibration_keywords` you pass, before `TREEDATA`;
- one tree record per row in the engine's default tree format, then `-999`, `PROCESS` and `STOP`.

Species are translated from FIA codes to the variant's alpha codes through the `FIAJSP` and `JSP` arrays in `config/calibrated/<variant>.json`, falling back to a short hand table when that file is absent. An FIA code with no mapping raises `KeyError`, so filter stems first; `filter_mapped_species(spcds, variant)` returns the indices to keep and a count of dropped codes. Crown ratio in percent becomes a class from 1 to 9 (`(cr + 9) // 10`, clipped). The tree number field is three digits wide, so keep each stand under 1,000 records.

The keyword file produced above ran through `fvs2py` in the documentation build and gave 21 summary rows from 2025 to 2125.

## Route 2: an SQLite input database

With the database route the keyword file names the database and the queries. The repository's benchmark template (`KEYFILE_TEMPLATE` in `calibration/python/bakuzis_100yr_comparison.py`) uses:

```text
STDIDENT
{stand_id}

DATABASE
DSNIN
{db_path}

DSNOUT
{db_path}

STANDSQL
SELECT *
FROM fvs_standinit
WHERE stand_id = '%StandID%'
ENDSQL

TREESQL
SELECT *
FROM fvs_treeinit
WHERE stand_id = '%StandID%'
ENDSQL

END
```

followed by output, cycle and calibration keywords and `PROCESS`/`STOP`. The tables are written with `pandas.DataFrame.to_sql`.

!!! warning "Unverified"
    The database route was not run for this page: the build used to write it lacked the SQLite C sources, including the column-name helper that decides how column names are matched. The repository writes lower-case table and column names (`fvs_standinit`, `stand_id`) while the Fortran readers compare against upper-case names (`STAND_ID`, `TREE_COUNT`); that combination is what the NE benchmark runner and microfvs use, but it was not checked here. The `fvs2py.FVS.add_trees` docstring also states that the in-process DBS read can leave example-tree arrays unpopulated and segfault; this was not reproduced.

### FVS_StandInit columns

The engine's reader (`src-converted/vdbsqlite/dbsstandin.f90`) recognises these columns, among others: `STAND_CN`, `STAND_ID`, `INV_YEAR`, `LATITUDE`, `LONGITUDE`, `REGION`, `FOREST`, `DISTRICT`, `COMPARTMENT`, `ECOREGION`, `LOCATION`, `HABITAT` or `PV_CODE`, `PV_REF_CODE`, `AGE`, `ASPECT`, `SLOPE`, `ELEVATION`, `ELEVFT`, `BASAL_AREA_FACTOR`, `INV_PLOT_SIZE`, `BRK_DBH`, `NUM_PLOTS`, `NONSTK_PLOTS`, `SAM_WT`, `STK_PCNT`, `DG_TRANS`, `DG_MEASURE`, `HTG_TRANS`, `HTG_MEASURE`, `MORT_MEASURE`, `SITE_SPECIES`, `SITE_INDEX`, `MAX_BA`, `MAX_SDI`, `MODEL_TYPE`, `PHYSIO_REGION`, `FOREST_TYPE`, `STATE`, `COUNTY`, plus fuel and photo-series columns.

The microfvs data model (`deployment/microfvs/microfvs/models.py`, `FvsStandInit`) makes six fields mandatory: `stand_id`, `variant`, `inv_year`, `basal_area_factor`, `inv_plot_size` and `brk_dbh`. The engine reader has no `VARIANT` case; microfvs uses that column to choose the executable.

The three design columns decide how `TREE_COUNT` is expanded. The FIA stand generator described below sets `basal_area_factor = -1.0`, `inv_plot_size = 1.0` and `brk_dbh = 999.0` so that every tree is on one fixed 1 acre plot and `TREE_COUNT` passes through as trees per acre. Its comments record the two errors this replaced: a basal area factor of 0 gave every tree at or above the break DBH zero expansion, and an `inv_plot_size` of 6.0 multiplied every record by six.

For Northeast stands, supply `ECOREGION`. NE switches the FIA volume, biomass and carbon routines on by default (`src-converted/ne/grinit.f90`, line 166), and without an ecodivision every run logs `FVS43 WARNING: FIAVBC REQUESTED BUT ECODIVISION CODE NOT RECOGNIZED OR MISSING` and falls back to `0000`.

### FVS_TreeInit columns

The engine's tree reader (`src-converted/dbsqlite/dbstreesin.f90`) recognises `PLOT_ID`, `TREE_ID`, `TREE_COUNT`, `HISTORY`, `SPECIES`, `DBH` or `DIAMETER`, `DG`, `HT`, `HTG`, `HTTOPK`, `CRRATIO`, `DAMAGE1` to `DAMAGE3`, `SEVERITY1` to `SEVERITY3`, `TREEVALUE`, `PRESCRIPTION`, `SLOPE`, `ASPECT`, `HABITAT` or `PV_CODE`, `TOPOCODE`, `SITEPREP`, `AGE`, `CULL`, `DECAYCD` and `WDLND_STEMS`. When `DBH` is zero or absent the reader uses `DIAMETER`.

| Column | Unit or coding |
|:---|:---|
| `TREE_COUNT` | Trees per acre represented by the record under the stand's design |
| `SPECIES` | FIA code or variant alpha code, read as text; a two digit numeric code is given a leading zero |
| `DIAMETER` / `DBH` | inches |
| `HT` | total height, feet |
| `CRRATIO` | live crown ratio, percent (microfvs validates 0 to 100) |
| `TREE_ID`, `PLOT_ID` | integers |

`TREE_ID` and `PLOT_ID` are read into default 4 byte integers (`dbstreesin.f90`, line 20), so an FIA control number (14 to 15 digits) cannot be used as a tree or plot identifier; `tools/fvs_input_guards.py` documents this.

The microfvs model (`FvsTreeInitRecord`) requires only `stand_id`, `species` and a `diameter` greater than zero.

!!! warning "Unknown species are remapped without an error"
    A species code that the variant does not carry is replaced by the variant's catch-all code, and the engine writes only a note to the main output (`NOTE: INPUT SPECIES CODE (...) WAS SET TO (...) FOR THIS PROJECTION.`, `src-converted/base/intree.f90`, line 421). `tools/fvs_input_guards.py` records that FIA codes 11, 20, 263 and 264, and the non-code 9999, all become `OH` in NE and give identical output. The run returns exit code 1 (notes only). Validate species codes before running.

## Building stands from FIA with fia_stand_generator

`calibration/python/fia_stand_generator.py` draws real FIA conditions into the stand and tree DataFrames used above. It was written for the Bakuzis benchmark and covers five variants through `VARIANT_STATES`: `ne` (CT, ME, MA, NH, NY, RI, VT), `acd` (ME, NH, VT), `pn` (OR, WA), `sn` (AL, FL, GA, MS, SC, TN) and `ie` (ID, MT, WA).

It reads per-state FIA CSV exports named `<ST>_COND.csv` and `<ST>_TREE.csv` from `fia_dir`, or from `FIA_DATA_DIR` (default `~/fia_data`) when `fia_dir` is not given.

```python
import sys
sys.path.insert(0, "calibration/python")

from fia_stand_generator import generate_real_stand
from inventory_keyfile import make_inventory_keyfile

cells = generate_real_stand("ne", "Medium", "Medium", n_plots=5, seed=42)
for i, (stand_df, tree_df, cond_row) in enumerate(cells):
    key = make_inventory_keyfile(stand_df, tree_df, f"S{i:03d}", "ne", num_cycles=20)
    open(f"runs/fia_{i:03d}.key", "w").write(key)
```

`generate_real_stand(variant, site_class, density_class, fia_dir=None, n_plots=5, seed=42)` returns up to `n_plots` tuples of `(stand_df, tree_df, cond_row)` and an empty list when no condition matches. It:

1. keeps COND rows with `COND_STATUS_CD == 1`;
2. keeps rows whose `SICOND` falls in the site bin (`SI_BINS`, feet: Low 30 to 50, Medium 50 to 70, High 70 to 100) and whose `BALIVE` falls in the density bin (`BA_BINS`, ft² ac⁻¹: Low 50 to 90, Medium 90 to 150, High 150 to 240), dropping rows with either value missing;
3. samples `n_plots` of them without replacement with `numpy.random.default_rng(seed)`;
4. loads live trees (`STATUSCD == 1`) for the sampled plots and, with `condition_trees(trees, cond_row)`, keeps only the trees of the sampled condition (`CONDID`) and divides each tree's `TPA_UNADJ` by the condition's share of that tree's own plot size (`MICRPROP_UNADJ` for microplot trees, `SUBPPROP_UNADJ` for subplot trees, `MACRPROP_UNADJ` for macroplot trees; `condition_expansion_divisor`), so that tree counts are per acre of condition, the basis on which FIA defines `BALIVE`. The plot size is read from `TPA_UNADJ` (74.965, 6.018 or 0.9997, within 1 percent); a tree that fits none of them, or whose proportion is missing, is divided by `CONDPROP_UNADJ` and counted in a log warning. Dividing every tree by `CONDPROP_UNADJ`, as the first fix did, left multi-condition catalog stands at 0.80 to 1.18 times `BALIVE` where a condition boundary crossed the microplot or macroplot differently from the subplot (defect D5).

`build_stand_init` fixes `inv_year` at 2000 for every stand so that trajectories share a calendar, sets latitude and longitude to 0, uses the design described above, and takes age, aspect, slope, elevation, site index, state, county and forest type from the condition row with defaults of 40 years, 0, 15 percent, 1,000 ft and site index 60 when values are missing. `build_tree_init` writes `tree_count` from `TPA_UNADJ` (rounded to 0.01), `species` from `SPCD`, `diameter` from `DIA` (missing becomes 5.0 in), `ht` from `HT` (missing becomes 40 ft) and `crratio` from `CR` (multiplied by 10 when the largest value is 9 or less, missing becomes 50, clipped to 15 to 99).

!!! warning "Identifiers"
    The generator derives each `stand_id` from the plot control number. Assign your own stand identifiers, as in the example, before sharing keyword files or outputs.

### Records with a null TPA_UNADJ

Live tree records whose `TPA_UNADJ` is null are dropped, with a logged warning giving the stand and the number of records (`build_tree_init`, lines 232 to 239). The code comments identify these as trees carried from a prior periodic inventory for growth, removal and mortality accounting (`P2A_GRM_FLG = Y`), which carry no area expansion and contribute nothing to `BALIVE`. Earlier versions filled them with 6.0 trees per acre; the code comments attribute 250 of 256 eastern seed-gate failures in the October 2026 catalog run to that fill, and `config/catalog/README.md` records that the catalog was regenerated without it with no change in any verdict.

### Condition proportion

!!! success "Fixed in PR 135 (October 30, 2026)"
    Before PR 135 the generator joined trees to the sampled condition by `PLT_CN` only and passed `TPA_UNADJ` through unchanged. On a plot with more than one condition, a condition therefore received the other conditions' trees, and a partial condition was seeded at `CONDPROP_UNADJ` times its density. The `condition_trees` docstring records that this affected 125 of the 360 NE and ACD stands in the Paper 1 ensemble (defect D3), and that the function mirrors the catalog runner `validation/bakuzis_stress_v2/catalog/conus_catalog_runs.py`. Keyword files generated before that date from multi-condition plots should be regenerated.

The repository's model-fitting dataset builders take a different approach and keep only conditions with `CONDPROP_UNADJ >= 0.75` (`calibration/R/30_build_conus_dataset.R`, line 16; `calibration/R/30b_build_htdbh_input.R`, lines 46 and 143).
