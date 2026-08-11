#!/bin/bash
#SBATCH --job-name=inspect
#SBATCH --account=PUOM0008
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect2_%j.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect2_%j.err

module load gcc/12.3.0
module load R/4.4.0
export R_LIBS_USER="/path/to/user/path"

Rscript -e "
library(brms)
cat(brms loaded OKn)

# DG model
dg <- readRDS(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/diameter_growth_fit.rds)
cat(n===== DG MODEL =====n)
cat(Class:, class(dg), n)

# Formula
cat(Formula:n)
print(dg[[formula]])

# Fixed effects
fe <- fixef(dg)
cat(nFixed effects:n)
print(fe)

# Random effects
re <- ranef(dg)
cat(nRandom effects groups:, paste(names(re), collapse=, ), n)
for (g in names(re)) {
  r <- re[[g]]
  cat(nGroup:, g, n)
  cat( Dims:, paste(dim(r), collapse= x ), n)
  levs <- rownames(r[,,1,drop=FALSE])
  cat( N levels:, length(levs), n)
  cat( First 15 levels:, paste(head(levs, 15), collapse=, ), n)
  cat( Dimnames[[3]]:, paste(dimnames(r)[[3]], collapse=, ), n)
  vals <- r[, Estimate, 1]
  cat( First 8 estimates:, paste(round(head(vals, 8), 4), collapse=, ), n)
  cat( Range:, round(min(vals), 4), to, round(max(vals), 4), n)
}

# Mortality model
cat(n===== MORTALITY MODEL =====n)
m <- readRDS(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/mortality_fit.rds)
cat(Formula:n)
print(m[[formula]])
fe_m <- fixef(m)
cat(nFixed effects:n)
print(fe_m)
re_m <- ranef(m)
for (g in names(re_m)) {
  r <- re_m[[g]]
  levs <- rownames(r[,,1,drop=FALSE])
  cat(nGroup:, g, n)
  cat( N levels:, length(levs), n)
  cat( First 15 levels:, paste(head(levs, 15), collapse=, ), n)
  vals <- r[, Estimate, 1]
  cat( First 8 estimates:, paste(round(head(vals, 8), 4), collapse=, ), n)
}

# HD model
cat(n===== HD MODEL =====n)
hd <- readRDS(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/height_diameter_fit.rds)
cat(Formula:n)
print(hd[[formula]])
fe_hd <- fixef(hd)
cat(nFixed effects:n)
print(fe_hd)
re_hd <- ranef(hd)
for (g in names(re_hd)) {
  r <- re_hd[[g]]
  levs <- rownames(r[,,1,drop=FALSE])
  cat(nGroup:, g, n)
  cat( N levels:, length(levs), n)
  cat( First 10 levels:, paste(head(levs, 10), collapse=, ), n)
  cat( Dimnames[[3]]:, paste(dimnames(r)[[3]], collapse=, ), n)
}

# Standardization
cat(n===== STANDARDIZATION PARAMS =====n)
std_file <- ${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/standardization_params.rds
if (file.exists(std_file)) {
  sp <- readRDS(std_file)
  cat(Class:, class(sp), n)
  if (is.data.frame(sp)) print(sp)
  if (is.list(sp) && !is.data.frame(sp)) {
    for (nm in names(sp)) {
      cat(nm, :, paste(sp[[nm]], collapse=, ), n)
    }
  }
} else {
  cat(NOT FOUNDn)
}

# Check training data columns and SPCD
cat(n===== DG TRAINING DATA =====n)
d <- read.csv(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/diameter_growth.csv, nrows=3)
cat(Columns:, paste(names(d), collapse=, ), n)

cat(n===== FILES IN NE OUTPUT =====n)
cat(paste(list.files(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/), collapse=n), n)

cat(n===== FILES IN NE DATA =====n)
cat(paste(list.files(${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/), collapse=n), n)
"
