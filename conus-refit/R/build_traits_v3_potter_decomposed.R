##=============================================================================
## build_traits_v3_potter_decomposed.R
##
## v3 trait file: decomposes Potter's vulnerability_score (composite) into
## the three orthogonal source components: climate_exposure (CE), sensitivity
## (S), low_adaptive_cap (LAC). These are nearly independent (CE-S corr =
## -0.01, CE-LAC corr = 0.08) so adding LAC and using CE in place of the
## composite gives the model strictly more information without redundancy.
##
## Also refines variety trait values for DF and LP using Potter framework:
##   - Coastal DF: lower CE (mild Pacific climate), slightly higher S
##   - Rocky Mountain DF: higher CE (continental extremes), lower S
##   - Shore pine: lower CE (coastal mild), higher S (narrow habitat)
##   - Rocky Mountain LP: higher CE (cold/fire-driven), moderate S
##
## New trait_cols (9 traits, was 8):
##   1. wood_specific_gravity
##   2. shade_tolerance_num
##   3. softwood
##   4. leaf_longevity_months
##   5. max_ht_m
##   6. max_dbh_cm
##   7. climate_exposure       <-- REPLACES vulnerability_score (composite)
##   8. sensitivity
##   9. low_adaptive_cap       <-- NEW (independent dimension)
##
## Resulting W matrix is 9 columns; existing 8-column v8/v5/v4 fits are
## NOT compatible. Drivers need trait_cols updated to use v3.
##=============================================================================

suppressPackageStartupMessages({ library(data.table) })

IN_FILE  <- "/users/PUOM0008/crsfaaron/fvs-conus/traits/species_traits_v2.rds"
OUT_FILE <- "/users/PUOM0008/crsfaaron/fvs-conus/traits/species_traits_v3.rds"

traits <- as.data.table(readRDS(IN_FILE))
cat(sprintf("Loaded v2 traits: %d species\n", nrow(traits)))

# --- Adjust variety rows using Potter framework intent ---------------------
# Original 202 lumped: vuln_score 53.49, CE 53.20, S 53.21, LAC 53.99 (assumed
# similar). Splitting by varieties' typical climate exposure:

# Coastal Douglas-fir: Pacific marine climate is climatically buffered
# (low CE), but the variety has narrow habitat specificity (higher S).
cat("\nAdjusting Coastal DF (2020) climate components:\n")
traits[SPCD == 2020L, `:=`(
  climate_exposure = 40.0,   # mild marine West coast, lower exposure
  sensitivity      = 60.0,   # specialized habitat, more sensitive
  low_adaptive_cap = 50.0    # similar to species-level LAC
)]
print(traits[SPCD == 2020L, .(SPCD, COMMON_NAME, climate_exposure, sensitivity, low_adaptive_cap)])

# Rocky Mountain DF: continental/drought-prone climate (high CE), but
# generalist across a wide range of mountain habitats (lower S).
cat("\nAdjusting Rocky Mountain DF (2021):\n")
traits[SPCD == 2021L, `:=`(
  climate_exposure = 65.0,   # continental extremes, droughts
  sensitivity      = 45.0,   # generalist mountain species, less sensitive
  low_adaptive_cap = 55.0    # somewhat constrained migration
)]
print(traits[SPCD == 2021L, .(SPCD, COMMON_NAME, climate_exposure, sensitivity, low_adaptive_cap)])

# Shore pine: coastal acid bogs/salt spray narrow habitat (high S),
# climate buffered but rising sea / coastal storms (moderate CE).
cat("\nAdjusting Shore pine (1080):\n")
traits[SPCD == 1080L, `:=`(
  climate_exposure = 45.0,
  sensitivity      = 55.0,
  low_adaptive_cap = 50.0
)]
print(traits[SPCD == 1080L, .(SPCD, COMMON_NAME, climate_exposure, sensitivity, low_adaptive_cap)])

# Rocky Mountain lodgepole: fire-dependent and cold-tolerant (high CE
# under warming), moderate S, established broad range (low LAC).
cat("\nAdjusting Rocky Mountain LP (1081):\n")
traits[SPCD == 1081L, `:=`(
  climate_exposure = 60.0,
  sensitivity      = 45.0,
  low_adaptive_cap = 45.0
)]
print(traits[SPCD == 1081L, .(SPCD, COMMON_NAME, climate_exposure, sensitivity, low_adaptive_cap)])

# --- Diagnostic check: report new VCC-component distribution -----------
cat("\n=== Potter component distribution (post-v3) ===\n")
for (col in c("climate_exposure","sensitivity","low_adaptive_cap")) {
  vals <- traits[[col]]
  vals <- vals[is.finite(vals)]
  cat(sprintf("  %-22s n=%d  median=%.1f  q05=%.1f  q95=%.1f\n",
              col, length(vals), median(vals),
              quantile(vals, 0.05), quantile(vals, 0.95)))
}

# --- Save -----------------------------------------------------------------
saveRDS(traits, OUT_FILE)
cat(sprintf("\nWrote: %s\n", OUT_FILE))
cat("\nDrivers need to be updated:\n")
cat("  trait_cols <- c('wood_specific_gravity', 'shade_tolerance_num',\n")
cat("                   'softwood', 'leaf_longevity_months',\n")
cat("                   'max_ht_m', 'max_dbh_cm',\n")
cat("                   'climate_exposure', 'sensitivity', 'low_adaptive_cap')\n")
cat("\nP_trait will increase from 8 to 9.\n")
