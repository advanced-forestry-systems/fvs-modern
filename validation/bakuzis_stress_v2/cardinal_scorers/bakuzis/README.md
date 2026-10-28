# Bakuzis matrix evaluation: four-arm CONUS realism benchmark

Benchmarks four model "arms" by the biological realism of their projected
stand-attribute relationships, scored against the OBSERVED Bakuzis matrix built
from FIA remeasurement pairs. Follows A. Weiskittel's existing Bakuzis
convention in `fvs-conus` (see "Convention" below).

All outputs land in `/fs/scratch/PUOM0008/crsfaaron/bakuzis/` (home quota is
full; scratch only). Run R with `module load gcc/12.3.0 R/4.4.0`.

## Convention (adopted, not invented)

Based on:
- `fvs-conus/BAKUZIS_UNCERTAINTY_README.md`
- `fvs-conus/R/16_bakuzis_matrix_figure.R`  (the stand-level state vars + panels)
- `fvs-conus/python/bakuzis_matrix_heatmap.py` (the small-multiples matrix grid)

Aaron's framework treats the Bakuzis matrix as the set of relationships among
stand-level state variables, with the structural (height-independent) variables
as the backbone and biomass as secondary. The state variables used here:

| Var | Meaning                         | Class       |
|-----|---------------------------------|-------------|
| BA  | stand basal area (m2/ha)        | structural  |
| QMD | quadratic mean diameter (cm)    | structural  |
| TPH | trees per hectare               | structural  |
| SDI | Reineke stand density index     | structural  |
| RD  | relative density (SDI/SDImax %) | structural  |
| AGB | aboveground biomass (Mg/ha)     | SECONDARY * |

\* AGB is down-weighted: the conus arms' AGB uses an NE-calibrated height
anchor, so it is reported but excluded from the headline realism index. Emphasis
is on BA, QMD, TPH, SDI, RD per the task brief.

### Matrix cells

Each cell is the log-log OLS slope of a canonical bivariate relationship (the
off-diagonal of the Bakuzis matrix), in the ecologically standard orientation:

| cell      | y ~ x     | reading                                   |
|-----------|-----------|-------------------------------------------|
| QMD_TPH   | QMD ~ TPH | self-thinning region                      |
| TPH_BA    | TPH ~ BA  | density-BA boundary                       |
| BA_QMD    | BA ~ QMD  | BA accrues with mean tree size            |
| SDI_QMD   | SDI ~ QMD | SDI vs mean size                          |
| SDI_BA    | SDI ~ BA  | SDI vs BA (near-1:1 in nature)            |
| BA_TPH    | BA ~ TPH  | BA vs density                             |
| RD_QMD    | RD ~ QMD  | relative density vs size                  |
| AGB_BA    | AGB ~ BA  | biomass-BA (secondary, ht-anchored)       |

Headline diagnostics (the classic Bakuzis "laws"):
- **Self-thinning boundary**: slope of the 99th-percentile upper edge of
  log(TPH) ~ log(QMD). Reineke reference -1.605.
- **SDI ceiling**: 99th-pct of SDI vs observed SDImax.

### Units

The observed remeasurement RDS is metric (BA m2/ha, QMD cm, TPH trees/ha, SDI,
RD %). The per-cycle arm CSVs are imperial; `score_arm.R` converts them to the
same metric units before computing any cell, so observed and arm cells are
directly comparable:
`BA *0.2295684`, `QMD *2.54`, `TPH` unchanged, `AGB *2.2417`,
`SDI = TPH*(QMD_cm/25)^1.605`, `RD = 100*SDI/SDImax`.

## Scoring metric (per cell)

```
slope_dev  = arm_slope - obs_slope          (signed)
sign_match = sign(arm_slope) == sign(obs_slope)   (right direction?)
cell_score = exp(-|slope_dev| / 0.5)        (1 = perfect, ->0 worse)
```
- **realism_index** = mean(cell_score) over the 7 structural cells (0..1).
- **sign_match_rate** = fraction of structural cells with the correct sign.
- **self_thin_score** = exp(-|arm_Reineke - obs_Reineke| / 0.3).
- AGB cell reported but excluded from realism_index.

## The four arms

| Arm                       | files glob                                                                                   | --config    |
|---------------------------|---------------------------------------------------------------------------------------------|-------------|
| engine default            | `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_engine_v2/conus_*_b*.csv`              | `default`   |
| engine calibrated         | `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_engine_v2/conus_*_b*.csv`              | `calibrated`|
| conus species-dependent   | `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq/conus_eq_*_b2_metrics.csv` | `conus_b2`  |
| conus species-free        | `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq/conus_eq_*_b1_metrics.csv` | `conus_b1`  |

## How to run

```bash
cd /fs/scratch/PUOM0008/crsfaaron/bakuzis
module load gcc/12.3.0 R/4.4.0

# 1. observed reference (run once; ~1 min on srun)
srun -A PUOM0008 -t 25 -c 4 --mem=16G Rscript build_observed_matrix.R

# 2. score each arm (run after that arm's jobs finish)
srun -A PUOM0008 -t 15 -c 2 --mem=8G Rscript score_arm.R \
  --arm-name engine_default \
  --glob '/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_engine_v2/conus_*_b*.csv' \
  --config default

srun -A PUOM0008 -t 15 -c 2 --mem=8G Rscript score_arm.R \
  --arm-name engine_calibrated \
  --glob '/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_engine_v2/conus_*_b*.csv' \
  --config calibrated

srun -A PUOM0008 -t 15 -c 2 --mem=8G Rscript score_arm.R \
  --arm-name conus_species_dependent_b2 \
  --glob '/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq/conus_eq_*_b2_metrics.csv' \
  --config conus_b2

srun -A PUOM0008 -t 15 -c 2 --mem=8G Rscript score_arm.R \
  --arm-name conus_species_free_b1 \
  --glob '/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq/conus_eq_*_b1_metrics.csv' \
  --config conus_b1
```

## Outputs

| file | content |
|------|---------|
| `observed_bakuzis_matrix.csv` | observed slope per (stratum x cell); stratum=ALL is the national target, plus EPA_L2_CODE strata |
| `observed_self_thinning.csv`  | observed Reineke boundary slope + SDImax |
| `observed_state_summary.csv`  | observed univariate distribution of each state var |
| `fig_observed_bakuzis_matrix.png` | observed log-log matrix figure with fitted slopes |
| `arm_<name>_bakuzis_cells.csv`    | per-cell arm slope, obs slope, deviation, sign_match, cell_score |
| `arm_<name>_bakuzis_summary.csv`  | one-row arm realism_index, sign_match_rate, self_thin_score, SDI ceiling |

To compare arms after all four are scored:
`cat arm_*_bakuzis_summary.csv` (or `rbind` in R / pandas) gives the headline
realism table across arms.

## Notes / caveats
- Observed strata use EPA_L2_CODE (present in the RDS). The arm CSVs are keyed by
  STAND_CN + VARIANT; the per-variant standinit carries a Bailey ECOREGION code,
  which is NOT the same code system as the RDS EPA codes. The scorer therefore
  benchmarks against the national (ALL) observed cells by default. Ecoregion-
  stratified scoring would require an EPA L1/L2 crosswalk for the standinit
  STAND_CN keys (follow-on).
- `score_arm.R` runs on whatever batches are present (partial arms OK), so it can
  be re-run as more batches land.
