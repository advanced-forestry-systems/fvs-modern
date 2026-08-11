#!/bin/bash
#SBATCH --job-name=inspect
#SBATCH --account=PUOM0008
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect_%j.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect_%j.err

module load gcc/12.3.0
module load R/4.4.0
export R_LIBS_USER="/path/to/user/path"

Rscript -e "
library(brms)

cat(\"\\n===== DG MODEL (NE) =====\\n\")
dg <- readRDS(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/diameter_growth_fit.rds\")
cat(\"Class:\", class(dg), \"\\n\")
cat(\"Formula:\\n\")
print(dg\$formula)
cat(\"\\nFixed effects:\\n\")
print(fixef(dg))
cat(\"\\nRandom effects groups:\", paste(names(ranef(dg)), collapse=\", \"), \"\\n\")
re <- ranef(dg)
for (g in names(re)) {
  cat(\"\\nGroup:\", g, \"\\n\")
  r <- re[[g]]
  cat(\"  Dims:\", paste(dim(r), collapse=\" x \"), \"\\n\")
  cat(\"  Level names (first 15):\", paste(head(rownames(r[,,1,drop=FALSE]), 15), collapse=\", \"), \"\\n\")
  cat(\"  N levels:\", nrow(r[,,1,drop=FALSE]), \"\\n\")
  cat(\"  Dimnames[[3]]:\", paste(dimnames(r)[[3]], collapse=\", \"), \"\\n\")
  cat(\"  First 5 Estimate values:\\n\")
  print(head(r[, \"Estimate\", , drop=FALSE], 5))
}

cat(\"\\n===== MORTALITY MODEL (NE) =====\\n\")
m <- readRDS(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/mortality_fit.rds\")
cat(\"Class:\", class(m), \"\\n\")
cat(\"Formula:\\n\")
print(m\$formula)
cat(\"\\nFixed effects:\\n\")
print(fixef(m))
cat(\"\\nRandom effects groups:\", paste(names(ranef(m)), collapse=\", \"), \"\\n\")
re_m <- ranef(m)
for (g in names(re_m)) {
  cat(\"\\nGroup:\", g, \"\\n\")
  r <- re_m[[g]]
  cat(\"  Level names (first 15):\", paste(head(rownames(r[,,1,drop=FALSE]), 15), collapse=\", \"), \"\\n\")
  cat(\"  N levels:\", nrow(r[,,1,drop=FALSE]), \"\\n\")
  cat(\"  First 10 Estimate values:\", head(r[, \"Estimate\", 1], 10), \"\\n\")
}

cat(\"\\n===== HD MODEL (NE) =====\\n\")
hd <- readRDS(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/height_diameter_fit.rds\")
cat(\"Class:\", class(hd), \"\\n\")
cat(\"Formula:\\n\")
print(hd\$formula)
cat(\"\\nFixed effects:\\n\")
print(fixef(hd))
cat(\"\\nRandom effects groups:\", paste(names(ranef(hd)), collapse=\", \"), \"\\n\")
re_hd <- ranef(hd)
for (g in names(re_hd)) {
  cat(\"\\nGroup:\", g, \"\\n\")
  r <- re_hd[[g]]
  cat(\"  Level names (first 10):\", paste(head(rownames(r[,,1,drop=FALSE]), 10), collapse=\", \"), \"\\n\")
  cat(\"  N levels:\", nrow(r[,,1,drop=FALSE]), \"\\n\")
  cat(\"  Dimnames[[3]]:\", paste(dimnames(r)[[3]], collapse=\", \"), \"\\n\")
}

cat(\"\\n===== STANDARDIZATION PARAMS =====\\n\")
std_file <- \"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/standardization_params.rds\"
if (file.exists(std_file)) {
  sp <- readRDS(std_file)
  cat(\"Class:\", class(sp), \"\\n\")
  print(sp)
} else {
  cat(\"File not found\\n\")
}

cat(\"\\n===== DG TRAINING DATA =====\\n\")
dg_data <- read.csv(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/diameter_growth.csv\", nrows=3)
cat(\"Columns:\", paste(names(dg_data), collapse=\", \"), \"\\n\")
dg_full <- read.csv(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/diameter_growth.csv\")
cat(\"N rows:\", nrow(dg_full), \"\\n\")
cat(\"N unique SPCD:\", length(unique(dg_full\$SPCD)), \"\\n\")
cat(\"SPCD values:\", paste(sort(unique(dg_full\$SPCD)), collapse=\", \"), \"\\n\")

cat(\"\\n===== MORTALITY TRAINING DATA =====\\n\")
m_data <- read.csv(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/mortality.csv\", nrows=3)
cat(\"Columns:\", paste(names(m_data), collapse=\", \"), \"\\n\")
m_full <- read.csv(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/mortality.csv\")
cat(\"N rows:\", nrow(m_full), \"\\n\")
cat(\"N unique SPCD:\", length(unique(m_full\$SPCD)), \"\\n\")

cat(\"\\n===== FILES IN NE OUTPUT =====\\n\")
print(list.files(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/NE/\"))

cat(\"\\n===== SLOPE/ELEV IN FIA COND TABLE =====\\n\")
fia_tree <- \"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/fia_tree.csv\"
if (file.exists(fia_tree)) {
  ft <- read.csv(fia_tree, nrows=3)
  cat(\"FIA tree columns:\", paste(names(ft), collapse=\", \"), \"\\n\")
} else {
  # Check what data files exist
  cat(\"FIA data files:\\n\")
  print(list.files(\"${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/NE/\"))
}
"
