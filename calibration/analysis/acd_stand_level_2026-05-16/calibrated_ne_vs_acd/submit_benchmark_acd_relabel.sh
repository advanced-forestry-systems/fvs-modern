#!/bin/bash
#SBATCH --job-name=fvs_bench_acd_relabel
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/fia_benchmark_acd_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/fia_benchmark_acd_%j.err

module purge
module load gcc/12.3.0 R/4.4.0 gdal/3.7.3 proj/9.2.1 geos/3.12.0

export FVS_PROJECT_ROOT="/users/PUOM0008/crsfaaron/fvs-modern-acdbridge"
export FIA_DATA_DIR="/fs/scratch/PUOM0008/crsfaaron/FIA"
export FVS_ACD_RELABEL=TRUE

cd ${FVS_PROJECT_ROOT}/calibration

Rscript R/19_fia_benchmark_engine.R 2>&1 | tee logs/fia_benchmark_acd_console.log
