#!/bin/bash
#SBATCH --job-name=cat_stress5 --account=PUOM0008 --time=01:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/stress5_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog
Rscript stress_catalog.R --long=catalog_long_v5.csv --traj=score_catalog_v5/catalog_20261029v5_scorecard_by_trajectory.csv --thresholds=fia_variant_thresholds_plot_acd.csv --gates=gates_v5/catalog_gates_by_arm.csv --out=stress_v5 --boot=500
echo STRESS_DONE
