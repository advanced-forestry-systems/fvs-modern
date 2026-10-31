#!/bin/bash
#SBATCH --job-name=cat_stress8 --account=PUOM0008 --time=01:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/stress8_%j.out
# Catalog v8 (D5 plot-size fix, PR 144): 144-cell rule grid and 500-draw stand bootstrap, same script and settings as stress_v6.slurm.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
C=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog; F=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030
cd $F
Rscript $C/stress_catalog.R --long=catalog_long_v8.csv --traj=score_catalog_v8/catalog_20261030v8_scorecard_by_trajectory.csv --thresholds=$C/fia_variant_thresholds_plot_acd.csv --gates=gates_v8/catalog_gates_by_arm.csv --out=stress_v8 --boot=500
echo STRESS_DONE
