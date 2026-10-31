# Outputs

A projection produces three kinds of output: the summary table that `fvs2py` returns as a DataFrame, the main output listing that FVS writes beside the keyword file, and any database tables that the keyword file requests. This page documents the summary table in detail.

## The summary table

`FVS.summary` reads the engine's per-cycle summary array (`IOSUM` in `src-converted/common/OUTCOM.f90`) through the `fvsSummary` API routine and returns it as a pandas DataFrame with these 20 columns, in this order (`deployment/fvs2py/fvs2py/constants.py`, `SUMMARY_COLS`). Values are integers: the engine rounds each per-acre quantity to the nearest whole unit before storing it (`src-converted/vbase/disply.f90`).

| Column | Content | Unit (engine native) |
|:---|:---|:---|
| `year` | Calendar year at the start of the cycle | year |
| `age` | Stand age | years |
| `tpa` | Live trees at the start of the cycle | trees ac⁻¹ |
| `tcuft` | Total cubic volume | ft³ ac⁻¹ |
| `mcuft` | Merchantable cubic volume | ft³ ac⁻¹ |
| `bdft` | Board foot volume | board ft ac⁻¹ |
| `rtpa` | Trees removed | trees ac⁻¹ |
| `rtcuft` | Total cubic volume removed | ft³ ac⁻¹ |
| `rmcuft` | Merchantable cubic volume removed | ft³ ac⁻¹ |
| `rbdft` | Board foot volume removed | board ft ac⁻¹ |
| `atba` | Basal area after treatment | ft² ac⁻¹ |
| `atccf` | Crown competition factor after treatment | dimensionless |
| `attopht` | Top (average dominant) height after treatment | ft |
| `prdlen` | Length of the period starting in this year (0 in the last row) | years |
| `acc` | Accretion over the period | ft³ ac⁻¹ yr⁻¹ |
| `mort` | Mortality over the period | ft³ ac⁻¹ yr⁻¹ |
| `sampwt` | Stand sampling weight, truncated to an integer; a weight below 1 is stored multiplied by 100,000 | as entered |
| `fortyp` | Forest type code | code |
| `sizecls` | Stand size class | code |
| `stkcls` | Stand stocking class | code |

The units follow the engine's own summary table header ("SUMMARY STATISTICS (PER ACRE OR STAND BASED ON TOTAL STAND AREA)", with `CU FT` and `BD FT` volume columns; `src-converted/vbase/sumout.f90`) and the variable definitions in `OUTCOM.f90` (accretion and mortality per acre per year). The removal columns are zero unless the keyword file schedules a harvest.

There is one row per cycle boundary, so a run of N cycles gives N + 1 rows; rows not yet simulated (all zero) are dropped. In a row where no trees were removed, `atba`, `atccf` and `attopht` hold the values at the start of the period, because the engine only distinguishes before and after treatment when there is a removal.

!!! success "Fixed in PR 143 (October 30, 2026)"
    The engine's summary array has 22 columns per cycle; the last two are sawlog cubic volume and sawlog cubic volume removed. `fvsSummary` copies all 22 (`src-converted/base/apisubs.f90`, `integer :: summary(22)`). Before this fix `fvs2py` allocated 20 integers per row, so it dropped the sawlog columns and wrote the last row's two extra integers past the end of its array. `FVS.summary` now returns all 22 columns: the 20 above, then `scuft` (sawlog cubic volume, ft³ ac⁻¹) and `rscuft` (sawlog cubic volume removed). It raises `RuntimeError` if the library reports a different row width. On 5 NE and 5 ACD stands the first 20 columns are identical before and after the fix (`validation/fvs2py_fixes_20261030/`: script `smoke_summary.py`, outputs `smoke_before.csv` and `smoke_after.csv`, comparison `smoke_comparison.csv`).

Volume columns depend on the volume equations compiled into the library. In the documentation build, which lacked the NVEL `.for` sources, all volume columns were zero; see the note on the [Quickstart](../getting-started/quickstart.md) page.

## Metric conversion

Neither `fvs2py` nor the summary table converts units. The library exports a `fvsUnitConversion` routine, and `fvs2py` binds it at load time, but no `FVS` method calls it. Convert after the run with the factors the engine itself uses (`src-converted/common/METRIC.f90`):

| From (engine) | To | Multiply by | Engine constant |
|:---|:---|:---|:---|
| ft² ac⁻¹ | m² ha⁻¹ | 0.2295643 | `FT2pACRtoM2pHA` |
| ft³ ac⁻¹ | m³ ha⁻¹ | 0.0699713 | `FT3pACRtoM3pHA` |
| trees ac⁻¹ | trees ha⁻¹ | 2.471 | `HAtoACR` |
| ft | m | 0.3048 | `FTtoM` |
| in | cm | 2.54 | `INtoCM` |

```python
import pandas as pd

def summary_to_metric(s: pd.DataFrame) -> pd.DataFrame:
    """Add SI columns to an fvs2py summary table."""
    out = s.copy()
    out["tph"] = s["tpa"] * 2.471                     # trees ha-1
    out["ba_m2ha"] = s["atba"] * 0.2295643            # m2 ha-1
    out["tvol_m3ha"] = s["tcuft"] * 0.0699713         # m3 ha-1
    out["mvol_m3ha"] = s["mcuft"] * 0.0699713         # m3 ha-1
    out["topht_m"] = s["attopht"] * 0.3048            # m
    out["acc_m3hayr"] = s["acc"] * 0.0699713          # m3 ha-1 yr-1
    out["mort_m3hayr"] = s["mort"] * 0.0699713        # m3 ha-1 yr-1
    return out
```

Because the engine rounds before `fvs2py` sees the values, converted figures carry the rounding of the imperial integers (for basal area, up to about 0.11 m² ha⁻¹). Board foot volume has no fixed metric equivalent and is left out.

## The main output listing

FVS writes a text listing named after the keyword file with the extension `.out` (`runs/quickstart.key` gives `runs/quickstart.out`). It contains the keyword echo, warnings and notes, per-cycle tree tables and the summary table with additional columns that `fvs2py` does not return: stand density index and quadratic mean diameter at the start of the period and after treatment, sawlog volumes and mean annual increment. Search it for `FVS` followed by two digits and `ERROR` or `WARNING` to find problems ([Troubleshooting](troubleshooting.md)). When GOMPIT is active the listing also records how many fitted species were matched (`GOMPMORT enabled: ...`).

## Database output

With an SQLite output database (`DATABASE` and `DSNOUT` in the keyword file), the `SUMMARY`, `TREELIDB` and `CARBREDB` keywords write the summary, tree list and carbon tables. microfvs reads such tables back (its `FvsOutputTableName` enumeration includes `fvs_atrtlist`, `fvs_calibstats` and `fvs_cases` among others).

!!! warning "Unverified"
    Database output was not produced for this page, because the build used to write it lacked the SQLite C sources; table and column names for database output are therefore not documented here.

## Per-tree state

`FVS.get_tree_attr` and `FVS.tree_table` are meant to read per-tree state in process but crash at commit fb4cf10; see [Running projections](running-projections.md#tree-lists).
