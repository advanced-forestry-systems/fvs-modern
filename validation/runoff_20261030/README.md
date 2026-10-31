# CONUS model run-off, October 30, 2026

This folder holds a pre-registered comparison of three model families, all run on the same 4,748 FIA catalog stands across 18 FVS variants:

- **A. FVS-modern engine**: default and calibrated parameters, each with native or GOMPIT mortality.
- **B. fvs-conus projector**: v2, plus v6 as a diagnostic only. Each runs b1 and b2 diameter growth with logit or GOMPIT mortality. Ingrowth is off.
- **C. Unmodified USDA FVS**: built from USDAForestService/ForestVegetationSimulator at commit b183c10f, using the upstream makefile.

Each contestant was scored two ways:

- **Long term:** 100 years at a 5-year step, scored with the Bakuzis stress test and the catalog gates.
- **Near term:** one cycle equal to each plot's own remeasurement interval, scored against observed time-2 survivor basal area, TPH, QMD and volume.

## Layout

| Path | Contents |
|---|---|
| `PREREG_runoff_20261030.md` | Pre-registration (sections 1 to 11) and amendments A1 to A3 (sections 12 to 14). |
| `runoff_report_20261030.md` / `.docx` / `.pdf` | Technical report (draft). |
| `figures/` | Fig a: near-term BA percent bias with 95% bootstrap intervals. Fig b: share of stands above FIA BA p99.9 by year 50, with Wilson intervals. Fig c: decision map. |
| `scripts/` | Run-off code, described below. |
| `slurm/` | The SLURM scripts as run on OSC Cardinal (account PUOM0008). |
| `tasks/` | Array task lists: 72 engine tasks and 288 projector tasks. |
| `results/` | Primary outputs (fixed calibrated code, amendment A3). |
| `results/sensitivity/` | Results from the pre-fix in-process calibrated path (Option 1, amendment A2). |

### Scripts

| Script | Role |
|---|---|
| `runoff_engine_runs.py` | Engine driver and projector input builder. |
| `runoff_nearterm_obs.R` | Builds the observed time-2 targets. |
| `runoff_b_project.sh` and `runoff_b_post.R` | Run the projector and post-process its output. |
| `make_zero_ingrowth.R` | Writes the zero-ingrowth lookup that turns ingrowth off. |
| `runoff_score_near.R` | Near-term scoring. |
| `runoff_long_prep.R` and `prep_fix.R` | Assemble the long-term tables. |
| `runoff_decide.R` | The pre-registered decision rule. |
| `d2_wrapper.py` | The audit of the in-process calibration SDI zero fill. |
| `runoff_figures.R` | Draws the figures. |

### Results

| File | Contents |
|---|---|
| `runoff_decision_by_variant.csv` | Decision table. |
| `catalog_gates_by_arm.csv` and `gate_table.csv` | Long-term gates. |
| `near_scores.csv` | Near-term bias, percent bias and RMSE with bootstrap CIs, for the common and all-stands sets. |
| `near_table_common.csv` | Compact near-term table. |
| `near_year0_check.csv` | Year-0 seeding check. |
| `fix_effect_calibrated_v2.csv` | Effect of the calibration fix. |
| `c_vs_a_default_native_full.csv` | USDA FVS against FVS-modern default native. |
| `d2_compare_by_variant.csv` | SDI zero-fill audit. |
| `bakuzis_scorecard_summary_*.csv` | Bakuzis relation summaries. |

No file here carries a plot identifier, a PLT_CN or coordinates. The per-stand joined files (`near_joined.csv`) and the run trajectories stay on Cardinal because they carry plot keys.

## How to rerun (OSC Cardinal)

Work root: `/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/`. The scripts expect this layout.

1. **Build the observations:**
   `Rscript scripts/runoff_nearterm_obs.R <pairs.rds> ../catalog/fia_variant_assignment.csv obs/nearterm_obs_frame.csv`
2. **Build USDA FVS:** use the upstream makefile with the FFLAGS and CFLAGS given in PREREG section 2C. Libraries go to `upstream_build/lib/`.
3. **Run the full set:** `bash slurm/submit_full.sh` submits the engines, projector inputs, projector array and scoring, chained with dependencies.
4. **Re-run calibrated near-term with the fixed FVS-modern code:**
   `sbatch slurm/engines_fix.slurm`, then `sbatch slurm/rescore_fix.slurm`. This is amendment A3 and needs `catalog/catalog_long_v6.csv`.
5. **Figures:** `Rscript scripts/runoff_figures.R`.

Inputs not in this folder:

- The FIA state tables (`$HOME/fia_data`).
- The CONUS remeasurement pairs RDS.
- The FVS-modern worktrees `wt-bakuzis-main` and `wt-fix-20261030`.
- The fvs-conus projector scripts in `fvs_stress/conus_eq_proj/`.

## Addendum, October 30, 2026: stand seeding and D5

The run-off seeded its stands with the catalog runner as it stood on October 30. That runner expanded each condition's trees by CONDPROP_UNADJ alone. FIA computes BALIVE with the condition's share of each plot size (MICRPROP_UNADJ, SUBPPROP_UNADJ, MACRPROP_UNADJ), and fvs-modern pull request 144 corrects the expansion to match (defect D5). Under the corrected expansion, 266 of the 4,748 catalog stands (5.6 percent) seed more than 2 percent differently. The catalog rerun on all four arms with it (catalog v8) changed no primary rule verdict, recommended arm or status. The run-off's conclusions are therefore expected to hold, but its per-stand numbers would shift slightly if the run-off were rerun with the corrected seeding. Memo Section 19 (`validation/bakuzis_stress_v2/20261027_fvs-modern_fvs-conus_bakuzis_assessment.md`) gives the details.
