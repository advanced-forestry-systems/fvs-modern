# Paper 1 Bakuzis ensemble, fourth regeneration (ens_full4, October 30, 2026)

Same design as ens_full3 (NE and ACD, 36 cells x 5 FIA stands, default + calibrated MAP + 100 posterior
draws, seed 42), run on Cardinal (array job 16014150, finish job 16014540) with the code in PR 135
(fix/draws-overlay-standgen-20261030): D1 draw keyword layout, D2 no MAP overlay on draws and no zero
SDImax, D3 condition-level tree join with CONDPROP_UNADJ expansion, D4 draw multipliers on the MAP basis.
Raw per-draw CSVs (about 771,000 rows) stay on Cardinal in bakuzis_v2_20261027/ens_full4/.

| path | content |
|---|---|
| agg/bakuzis_uncertainty_summary_long.csv | legacy summary; posterior_band rows pool draws x stands (band_basis column) |
| agg/bakuzis_uncertainty_stand_bands.csv | per stand: 2.5/50/97.5 percent of its own 100 draws, MAP, default (parameter-only band) |
| agg/bakuzis_uncertainty_cell_bands.csv | per cell: median within-stand band, and between-stand spread of the 5 MAP values |
| agg/bakuzis_benchmark_wide.csv, agg/bakuzis_laws_compliance.csv | as before |
| ens_full4_seedgate.csv | year-0 BA of every stand vs its BALIVE selection bin |
| score_ens_full4/, score_ens_full4b/ | stress scorer, variant and species-group thresholds |
| ens4_numbers.R, .out, ens4_numbers_old_vs_new.csv | every Section 4.6 number, ens_full3 next to ens_full4 |
| figures/ | Paper 1 figures; trajectory ribbon is the within-stand band, band decomposition shows both |
| confirm/ | real-stand checks of the fixes (keyword echo, MAP overlay, BALIVE reproduction, draw medians vs MAP) |
| scripts/ | slurm and R scripts, job logs |
| md5sums.txt | md5 of every file here |

No plot identifiers or coordinates are in any file.
