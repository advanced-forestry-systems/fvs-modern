#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J ne_ab_sdicon
#SBATCH -t 03:00:00
#SBATCH -N 1 -n 8 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/ne_ab_%j.log
module load gcc/12.3.0 R/4.4.0
P=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
echo '=== BASELINE (--noconstrain) ==='
Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=300 --seed=7 --horizon=100 --cyclelen=5 --noconstrain --outdir=$P/ne_ab_baseline 2>&1 | grep -E 'CONSTRAINT|stands projected|year 0|PROJ_YEAR|Wrote|Error|error' 
echo '=== CONSTRAINED (default stand-level) ==='
Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=300 --seed=7 --horizon=100 --cyclelen=5 --outdir=$P/ne_ab_constrained 2>&1 | grep -E 'CONSTRAINT|stands projected|year 0|PROJ_YEAR|Wrote|Error|error'
echo NE_AB_DONE
