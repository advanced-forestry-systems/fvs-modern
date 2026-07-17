#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v3_smoke
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_smoke/v3_smoke_%j.log
module load gcc/12.3.0 R/4.4.0
P=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
O=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_smoke
echo "===== A: constrained, annual (cyclelen=1) ====="
Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=200 --cyclelen=1 --outdir=$O/con_annual
echo "===== B: NO constraint, annual (v2 ad hoc ramp) ====="
Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=200 --cyclelen=1 --noconstrain --outdir=$O/noncon_annual
echo "===== C: constrained, 5yr step (path-invariance vs A) ====="
Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=200 --cyclelen=5 --outdir=$O/con_5yr
echo ALL_DONE
