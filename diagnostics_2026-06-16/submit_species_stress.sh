#!/bin/bash
#SBATCH --job-name=species_stress
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/species_stress_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/species_stress_%j.err
module load gcc/12.3.0
module load gdal/3.7.3 geos/3.12.0 proj/9.2.1
module load R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
Rscript --vanilla species_stress.R
