#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 24:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/v3arr_%A_%a.out
## Constrained, annualized v3 projection: 19 variants x 2 modes = 38 tasks.
## Idempotent per (variant,mode); stand-level density + BA constraint ON; fitted coefficients.
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/out_conus_eq_v3; MAN=$D/eq_manifest.tsv
TASK=${SLURM_ARRAY_TASK_ID}
LINE=$(awk -v t="$TASK" '$1==t{print;exit}' "$MAN")
VAR=$(echo "$LINE"|cut -f2); MODE=$(echo "$LINE"|cut -f3)
CFG=$([ "$MODE" = "dependent" ] && echo conus_b2 || echo conus_b1)
MET=$OUT/conus_eq_$(echo $VAR|tr 'A-Z' 'a-z')_${CFG}_metrics.csv
mkdir -p "$OUT"
if [ -s "$MET" ] && [ "$(wc -l <"$MET")" -gt 1 ]; then echo "already done; skip"; exit 0; fi
module load gcc/12.3.0 R/4.4.0
if [ "$MODE" = "free" ] && [ ! -f "$OUT/b1_re_means.rds" ]; then
  for src in "$D/b1_re_means.rds" "$D/out_conus_eq/b1_re_means.rds"; do [ -f "$src" ] && cp "$src" "$OUT/b1_re_means.rds" && break; done
fi
Rscript --vanilla $D/conus_eq_projector_v3.R --variant=$VAR --mode=$MODE --nstands=0 \
  --cyclelen=1 --st_rate=0.4525 --st_exp=2.73 --ba_g1=0.036 --ba_k=0.0529 --outdir="$OUT"
echo "TASK_DONE $VAR $MODE"
