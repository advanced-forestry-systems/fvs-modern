#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 04:00:00
#SBATCH -N 1 -n 4
#SBATCH --mem=100G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_agb_%A_%a.out
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
MAN=$D/greg_agb_manifest.tsv
TASK=${SLURM_ARRAY_TASK_ID}
VAR=$(awk -v t="$TASK" '$1==t{print $2}' "$MAN")
MET=$D/out_conus_eq_greg/conus_eq_${VAR}_greg_metrics.csv
TRE=$D/out_conus_eq_greg/conus_eq_${VAR}_greg_treelists.csv
echo "task $TASK: variant=$VAR"
module load gcc/12.3.0 python/3.12
export PYTHONNOUSERSITE=1
export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python $D/conus_eq_agb.py --metrics "$MET" --treelists "$TRE"
echo "GREG_AGB_DONE $VAR"
