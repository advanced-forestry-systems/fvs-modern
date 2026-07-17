#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 12:00:00
#SBATCH -N 1 --ntasks=1 --cpus-per-task=48
#SBATCH --mem=384G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_sn_par_%j.out
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=$D/out_conus_eq_greg; MET=$OUT/conus_eq_sn_greg_metrics.csv
mkdir -p "$OUT"
if [ -s "$MET" ] && [ "$(wc -l <"$MET")" -gt 1 ]; then echo "already done SN; skip"; exit 0; fi
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla $D/conus_eq_projector_greg_par.R --variant=SN --nstands=0 --outdir="$OUT"
echo SN_PAR_DONE
