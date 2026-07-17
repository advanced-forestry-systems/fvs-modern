#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 06:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/eqb2gompit_%A_%a.out
## CONUS+AK fvs-conus equation projection, GOMPIT mortality arm (arm 5).
## 19 variants x 1 mode (dependent DG = conus_b2, isolates Greg's mortality).
## Mirrors conus_eq_array.sbatch; output -> out_conus_eq_b2_gompit/, CONFIG=conus_b2_gompit.
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/out_conus_eq_b2_gompit
MAN=$D/eq_gompit_manifest.tsv
TASK=${SLURM_ARRAY_TASK_ID}
LINE=$(awk -v t="$TASK" '$1==t{print; exit}' "$MAN")
VAR=$(echo "$LINE" | cut -f2)
MODE=$(echo "$LINE" | cut -f3)
CFG=conus_b2_gompit
MET=$OUT/conus_eq_$(echo $VAR | tr 'A-Z' 'a-z')_${CFG}_metrics.csv
echo "task $TASK: variant=$VAR mode=$MODE cfg=$CFG (gompit)"
mkdir -p "$OUT"
## idempotency: skip if metrics already exist with >1 row and AGB filled
if [ -s "$MET" ] && [ "$(wc -l < "$MET")" -gt 1 ] && head -2 "$MET" | tail -1 | grep -qE ',[0-9.]+'; then
  echo "already done; skipping"; exit 0
fi
module load gcc/12.3.0 R/4.4.0
Rscript "$D/conus_eq_projector_v2.R" --variant=$VAR --mode=dependent --mort=gompit \
    --cch_module="$D/cch_module.R" --nstands=0 --outdir="$OUT"
module load python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
TRE=$OUT/conus_eq_$(echo $VAR | tr 'A-Z' 'a-z')_${CFG}_treelists.csv
python "$D/conus_eq_agb.py" --metrics "$MET" --treelists "$TRE"
echo "TASK_DONE $VAR gompit"
