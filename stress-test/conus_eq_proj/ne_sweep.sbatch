#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J ne_sweep
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 8 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/ne_sweep_%j.log
module load gcc/12.3.0 R/4.4.0
P=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
run(){ tag=$1; shift; od=$P/sweep_$tag; mkdir -p $od; Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=300 --seed=7 --horizon=100 --cyclelen=5 --outdir=$od "$@" >/dev/null 2>&1; echo "$tag -> $od"; }
run default
run gentle  --st_rate=0.15 --st_exp=3.5
run gentler --st_rate=0.08 --st_exp=4.0
run baG1hi  --st_rate=0.15 --st_exp=3.5 --ba_g1=0.045
echo SWEEP_DONE
