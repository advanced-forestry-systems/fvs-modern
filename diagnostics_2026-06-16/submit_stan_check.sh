#!/bin/bash
#SBATCH --job-name=stan_check
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/stan_check_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/stan_check_%j.err
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
cd ~/fvs-conus
Rscript -e 'library(cmdstanr); m<-cmdstan_model("stan/dg_kuehne2022_speciesfree_v8_forest_eco.stan", compile=TRUE); cat("COMPILE_OK\n"); cat("variables:\n"); print(m$variables()$parameters[c("z_FT_raw","sigma_FT")])' 2>&1 | tail -20
