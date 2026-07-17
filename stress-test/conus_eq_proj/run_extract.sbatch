#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J emt_td
#SBATCH -t 00:30:00
#SBATCH -N 1 -n 4 --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/emt_td_%j.out
module load gcc/12.3.0 R/4.4.0
SP=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids
export LD_LIBRARY_PATH=$SP/proj/gcc/12.3.0/9.2.1-iqydouj/lib64:$SP/gdal/gcc/12.3.0/3.7.3-jfsjqp3/lib64:$SP/geos/gcc/12.3.0/3.12.0-crfkwd5/lib64:$LD_LIBRARY_PATH
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript extract_emt_td.R
