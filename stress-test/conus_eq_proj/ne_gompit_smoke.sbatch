#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 4
#SBATCH --mem=24G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/ne_gompit_smoke_%j.out
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/ne_gompit_smoke
mkdir -p "$OUT" "$D/logs"
module load gcc/12.3.0 R/4.4.0
echo "===== GOMPIT arm (NE, 40 stands, seed 7) ====="
Rscript "$D/conus_eq_projector_v2.R" --variant=NE --mode=dependent --mort=gompit --nstands=40 --seed=7 --outdir="$OUT"
echo "===== LOGIT arm (NE, 40 stands, seed 7) [baseline] ====="
Rscript "$D/conus_eq_projector_v2.R" --variant=NE --mode=dependent --mort=logit  --nstands=40 --seed=7 --outdir="$OUT"
echo "ALL_SMOKE_DONE"
