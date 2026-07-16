#!/bin/bash
#SBATCH --job-name=fvs_overstory
#SBATCH --account=PUOM0008
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/logs/overstory_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/logs/overstory_%j.err
module purge
module load gcc/12.3.0 R/4.4.0 gdal/3.7.3 proj/9.2.1 geos/3.12.0
export FVS_PROJECT_ROOT="/users/PUOM0008/crsfaaron/fvs-modern"
export FVS_FIA_DATA_DIR="/users/PUOM0008/crsfaaron/FIA"
export FVS_ACD_RELABEL=TRUE
export FVS_ACD_FOOTPRINT_STATES="23,33,50"
cd "$FVS_PROJECT_ROOT/calibration"
Rscript R/19_fia_benchmark_engine_overstory.R 2>&1 | tee logs/overstory_console.log
