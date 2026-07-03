#!/bin/bash
#SBATCH --job-name=fourway_agg
#SBATCH --time=00:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourway_agg_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourway_agg_%j.err
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
Rscript --vanilla calib_4way_aggregate.R 2>&1 | tail -40
