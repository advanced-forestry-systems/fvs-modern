#!/bin/bash
#SBATCH --job-name=ro_full_b --account=PUOM0008 --time=06:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=12G
#SBATCH --array=0-287%60
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/logs/b_%A_%a.out
# Run-off full run, contestant B: fvs-conus projector v2 and v6 (ingrowth off), 18 variants x 4 arms x {long, near}.
# Task list full/tasks_b.tsv: idx VAR ARM MODE VER. Arms are labelled conus_v2_* / conus_v6_* by runoff_b_post.R.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full
read I V A M VER <<< "$(awk -v t=$SLURM_ARRAY_TASK_ID '$1==t' $F/tasks_b.tsv)"
echo ">>> B task $I $VER $V $A $M $(date)"
bash $RO/scripts/runoff_b_project.sh $V $A $M $F/b_inputs/$M $F/b_out $VER
echo "rc=$? $(date)"
