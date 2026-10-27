#!/bin/bash
#SBATCH --job-name=bakuzis6
#SBATCH --account=PUOM0008
#SBATCH --time=02:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis/bakuzis6_%j.log
module load gcc/12.3.0 R/4.4.0
B=/fs/scratch/PUOM0008/crsfaaron
cd $B/bakuzis

# Arm: engine default
Rscript score_arm.R --arm-name engine_default \
  --glob "$B/fvs_stress/out_conus_engine_v4/conus_*_b*.csv" \
  --config default --out $B/bakuzis

# Arm: engine calibrated
Rscript score_arm.R --arm-name engine_calibrated \
  --glob "$B/fvs_stress/out_conus_engine_v4/conus_*_b*.csv" \
  --config calibrated --out $B/bakuzis

# Arm: greg faithful
Rscript score_arm.R --arm-name greg_faithful \
  --glob "$B/fvs_stress/conus_eq_proj/out_conus_eq_greg/conus_eq_*_greg_metrics.csv" \
  --config greg --out $B/bakuzis

echo "BAKUZIS6 DONE"
