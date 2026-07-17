#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/b1gompit_smoke_%j.out
# b1 (species-free) growth + gompit mortality, NE smoke. Mirrors conus_eq_gompit_array.sbatch
# but with --mode=free (conus_b1) instead of dependent (conus_b2). Completes the b1 x gompit
# cell that the six-arm test never ran (only b2 x gompit = the "gompit" arm exists).
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/out_conus_eq_b1_gompit_smoke
mkdir -p "$OUT" "$D/logs"
module load gcc/12.3.0 R/4.4.0
Rscript "$D/conus_eq_projector_v2.R" --variant=NE --mode=free --mort=gompit \
    --cch_module="$D/cch_module.R" --nstands=1500 --outdir="$OUT"
module load python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
MET=$OUT/conus_eq_ne_conus_b1_gompit_metrics.csv
TRE=$OUT/conus_eq_ne_conus_b1_gompit_treelists.csv
python "$D/conus_eq_agb.py" --metrics "$MET" --treelists "$TRE"
echo "SMOKE_DONE; rows=$(wc -l < "$MET")"
