#!/bin/bash
#SBATCH --job-name=ingrowth_check
#SBATCH --time=00:12:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/ingrowth_check_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/ingrowth_check_%j.err
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
cd ~/fvs-conus
Rscript -e 'library(cmdstanr); m<-cmdstan_model("stan/ingrowth_negbinom_v3_forest_eco.stan", compile=TRUE); cat("INGROWTH_COMPILE_OK\n"); print(names(m$variables()$parameters)[grepl("FT", names(m$variables()$parameters))])' 2>&1 | tail -8
