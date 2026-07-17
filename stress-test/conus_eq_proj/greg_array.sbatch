#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 12:00:00
#SBATCH -N 1 -n 8
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_arr_%A_%a.out
## Faithful Greg arm (CONFIG=greg): 19 variants. Idempotent per variant.
## Memory is set at submit time via --mem on the sbatch CLI per tier (SN 384G; LS/NE/CR/CS 248G; rest 128G).
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/out_conus_eq_greg; MAN=$D/greg_manifest.tsv
TASK=${SLURM_ARRAY_TASK_ID}
VAR=$(awk -v t="$TASK" '$1==t{print $2;exit}' "$MAN")
MET=$OUT/conus_eq_$(echo $VAR|tr 'A-Z' 'a-z')_greg_metrics.csv
mkdir -p "$OUT"
if [ -s "$MET" ] && [ "$(wc -l <"$MET")" -gt 1 ]; then echo "already done $VAR; skip"; exit 0; fi
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla $D/conus_eq_projector_greg.R --variant=$VAR --nstands=0 --outdir="$OUT"
echo "GREG_TASK_DONE $VAR"
