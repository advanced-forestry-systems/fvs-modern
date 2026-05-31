# SILC strata regenerated on the 5-type x 2-density break

SILC requested a revised stratification for the multi-model SILC
benchmark. Replaces the prior 11 byStrata stand grouping with a
landscape grouping of 10 cells:

|  | A+B (high density) | C+D (low density) |
| --- | --- | --- |
| Cedar               | 28   plots | 168  plots |
| Hardwood            | 534  plots | 1063 plots |
| Mixedwood (HS + SH) | 348  plots | 510  plots |
| Commercial Softwood | 662  plots | 1018 plots |
| Other Softwood      | 35   plots | 121  plots |

SILC stand-ID coding (e.g. "S3B-N"):

* Alphabetic prefix: forest type
  * `C`  -> Cedar
  * `H`  -> Hardwood
  * `HS` -> Hardwood-leading Mixedwood   } folded to Mixedwood
  * `SH` -> Softwood-leading Mixedwood   }
  * `S`  -> Commercial Softwood
  * `OS` -> Other Softwood
* 3rd letter A/B/C/D: density class. A+B = high, C+D = low.

## Files in this directory

| File | What it is | Status |
| --- | --- | --- |
| `silc_strata_5x2_regen.R`           | regenerator script | local |
| `silc_strata_5x2_mapping.csv`       | 79 matrix-stand strata mapped to (forest_type, density_class) | local |
| `silc_strata_5x2_rollup.csv`        | plot counts rolled up to the 10 cells | local |
| `silc_strata_5x2_plot_counts.png`   | plot-count figure (5x2 bars) | local |
| `silc_strata_5x2_AGM_trajectories.csv` | AGM strata-mean trajectories on the new break | local |
| `silc_strata_5x2_AGM_BA.png`        | BA 100-yr trajectory by 5x2 cell, AGM only | local |
| `silc_strata_5x2_AGM_TPA.png`       | TPA 100-yr trajectory | local |
| `silc_strata_5x2_AGM_QMD.png`       | QMD 100-yr trajectory | local |
| `silc_strata_5x2_AGM_NetCords.png`  | net merchantable cords 100-yr trajectory | local |

## Data status by model

### AGM (AcadianGY) — local, complete

Source: `silc_extracted/GrownDB_byStrata_ALL.csv` (144 MB, 21 cycles
from year 2023 to 2123, 11 SILC byStrata stands). All 11 stands cover
6 of the 10 new cells:

| Cell | n stands | byStrata IDs |
| --- | --- | --- |
| Cedar / C+D                 | 1 | C4C-N |
| Hardwood / A+B              | 1 | H3B-N |
| Hardwood / C+D              | 3 | H3C-N, H3D-N, H4C-N |
| Mixedwood / C+D             | 1 | HS3C-N |
| Commercial Softwood / A+B   | 2 | S2A-N, S3B-N |
| Commercial Softwood / C+D   | 3 | S3C-N, S3D-N, S4C-N |

Year-0 cross-check (DBH >= 4.6 in) reproduces SILC mature-Acadian
range: BA 69-164 ft^2/ac, TPA 168-394 trees/ac, QMD 6.8-9.5 in.

Empty cells (no AGM byStrata stand): Cedar / A+B, Mixedwood / A+B,
Other Softwood / A+B, Other Softwood / C+D. These appear as
"no AGM data" panels in the figures.

### OSM-ACD — NOT in current local checkout

The local `OSM.TreeListProjections.csv` upload (135 MB) is the
**FIA paired-plot** output (2437 numeric SurveyIDs, year 2014-2019,
cycle 0-5). That is the FIA accuracy benchmark dataset, not the SILC
100-year projection. Confirmed via column inspection:

* `SurveyID` is purely numeric (FIA plot CN), no SILC matrix-stand IDs
  like `S3B-N` are present
* `Year` 2014-2019 only
* `Cycle` 0-5 (5 year FIA remeasurement)

The OSM-ACD 100-year SILC trajectories were generated in prior
sessions on Cardinal and were aggregated into the `silc_v25..v28`
analysis tree that does not mount in this checkout.

### FVS-NE / FVS-ACD — NOT in current local checkout

Same situation as OSM-ACD: prior-session Cardinal output, aggregated
locally into the missing `silc_v25..v28` tree. None of the four-model
SILC trajectory CSVs are present in `repos--fvs-modern/`.

## Year-100 outcomes and reconciliation diagnostics

* `silc_strata_5x2_year100_outcomes.csv` -- per-cell year-0 and
  year-100 BA, TPA, QMD, NetCords with 100-year growth factors
* `silc_strata_5x2_year100_growth.png` -- 5x2 grouped bar of BA and
  merch-cords growth factors by cell, n/a cells flagged in red
* `silc_strata_5x2_reconciliation_table.csv` -- the 11 byStrata
  stands -> 10-cell crosswalk
* `silc_strata_5x2_reconciliation.png` -- one-page Sankey-like
  diagram showing the prior 11-stand grouping flowing into the new
  10-cell break with empty cells flagged

Headline year-100 BA growth factors (AGM only):

| Cell | n stands | BA factor | NetCords factor |
| --- | --- | --- | --- |
| Cedar / C+D                | 1 | 1.28 | 1.66 |
| Hardwood / A+B             | 1 | 2.15 | 2.62 |
| Hardwood / C+D             | 3 | 2.66 | 3.41 |
| Mixedwood / C+D            | 1 | 2.17 | 2.85 |
| Commercial Softwood / A+B  | 2 | 1.61 | 2.25 |
| Commercial Softwood / C+D  | 3 | 1.93 | 2.54 |

The C+D (low-density) cells project larger growth factors than
A+B because they start from a smaller base and have more room
to fill in to maturity. Cedar has the lowest growth factor
(1.28x) reflecting slow growth of mature C4C-N white-cedar stand.

## Relative-change and uniform-axis figure variants

* `silc_strata_5x2_AGM_relBA.png` -- BA trajectories normalized to
  year-0 value per cell; y in growth-factor units, horizontal at 1.0
* `silc_strata_5x2_AGM_relNetCords.png` -- merch cords growth factor
  trajectories
* `silc_strata_5x2_AGM_BA_unifrow.png` -- BA in absolute units with
  y-axis fixed within each density row so cross-forest-type
  comparison is readable within the row
* `silc_strata_5x2_AGM_NetCords_unifrow.png` -- same for net cords

## CFI empirical evaluation

`silc_cfi/` (see its own `README.md`) anchors the strata work with
actual Seven-Islands remeasurement data: 10 fixed-plot CFI plots
measured 1981-2000 (24 reliable remeasurement intervals).

AcadianGY 12.3.9 on the 17 routine-growth CFI pairs:
* BA bias: +5.96%
* RMSE: 7.48 ft^2/ac
* R^2: 0.91
* Naive FIA-prior baseline: +3.94% bias, RMSE 7.03

The CFI evaluation maps cleanly to the 5x2 break:

| Cell | CFI plots |
| --- | --- |
| Cedar / A+B           | 1106 |
| Hardwood / A+B        | 1107 |
| Hardwood / C+D        | 1104 |
| Mixedwood / A+B       | 1101, 1102, 1109 |
| Mixedwood / C+D       | 1100 |
| Commercial SW / C+D   | 1105, 1108 |

That gives CFI coverage in 6 of the 10 cells, including the
Cedar / A+B cell which AGM doesn't cover. The two coverages
(AGM trajectories at long horizon; CFI empirical at short horizon)
are complementary.

## To extend to all four models on the new break

Two paths:

1. Restore the prior `silc_v25..v28` analysis tree. The aggregated
   per-stand-per-year-per-model CSVs there can be remapped onto the
   new 5-type / 2-density key using `silc_strata_5x2_mapping.csv`
   without re-running anything. The regenerator script already
   contains the mapping function `parse_stand_id()`.

2. Re-run the four-model SILC projection on Cardinal and pull the
   aggregated CSVs down. The submit scripts existed in
   `calibration/osc/` (`submit_silc_*.sh`); the SILC tree carried
   the aggregator (`R/silc_aggregate_by_stratum.R`).

Once any per-model per-stand-per-year CSV is available locally, the
five-column join `c(model, stand_id, year, ba_ft2ac, ...)` merged
against `silc_strata_5x2_mapping.csv` reproduces the 10-cell rollup
in seconds; the figure code in `silc_strata_5x2_regen.R` just needs
the loop extended over `unique(model)`.

## Year-0 reasonableness check (AGM only)

```
StandID   TPA     BA   QMD       Cell
C4C-N   359.8  161.6  9.08     Cedar / C+D
H3B-N   234.3   97.2  8.72     Hardwood / A+B
H3C-N   225.7   75.9  7.85     Hardwood / C+D
H3D-N   229.1   68.7  7.42     Hardwood / C+D
H4C-N   167.7   79.2  9.30     Hardwood / C+D
HS3C-N  267.0   97.7  8.19     Mixedwood / C+D
S2A-N   347.0   88.4  6.83     Commercial Softwood / A+B
S3B-N   394.3  164.3  8.74     Commercial Softwood / A+B
S3C-N   263.4  100.6  8.37     Commercial Softwood / C+D
S3D-N   328.5   89.7  7.07     Commercial Softwood / C+D
S4C-N   283.3  138.1  9.45     Commercial Softwood / C+D
```

Filter: DBH >= 4.6 in (the BRK_DBH used in SILC StandInit).
All eleven stands fall in the mature-Acadian Matrix envelope.
