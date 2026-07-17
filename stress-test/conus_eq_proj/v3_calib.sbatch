#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v3_calib
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_calib/v3_calib_%j.log
module load gcc/12.3.0 R/4.4.0
P=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
O=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_calib
for R in 0.30 0.38 0.4525; do
  echo "===== ST_RATE=$R (ST_EXP=2.73) ====="
  Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=300 --cyclelen=1 --st_rate=$R --outdir=$O/rate_$R 2>/dev/null | grep -E "Realized self-thinning|PROJ_YEAR|1:|3:"
done
echo CALIB_DONE
