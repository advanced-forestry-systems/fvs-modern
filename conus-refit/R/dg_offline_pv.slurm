#!/bin/bash
#SBATCH --job-name=dg_offline_pv
#SBATCH --account=PUOM0008
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-conus/logs/dg_offline_pv_%j.out
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load gdal/3.7.3 geos/3.12.0 proj/9.2.1
module load R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-conus
Rscript R/dg_default_recal_per_variant.R
