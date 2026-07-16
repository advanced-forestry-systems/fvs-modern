##=============================================================================
## build_traits_v2_df_split.R
##
## Adds variety/subspecies split rows to species_traits.rds for species where
## a single SPCD lumps biologically and operationally distinct varieties:
##
##   1. Douglas-fir (SPCD 202):
##      - 2020 = Coastal (var. menziesii)   — EPA L1 = 7
##      - 2021 = Rocky Mountain (var. glauca) — else
##
##   2. Lodgepole pine (SPCD 108):
##      - 1080 = Shore pine (var. contorta)   — EPA L1 = 7
##      - 1081 = Rocky Mountain/Sierra lodgepole (var. latifolia + murrayana) — else
##
## Trait values informed by USDA Silvics Manual (Burns & Honkala 1990).
## Working max_ht_m and max_dbh_cm use typical mature working-forest values,
## not champion-tree extremes. Other traits adjusted modestly to reflect
## variety-level biology where defensible.
##
## The original SPCD 202 / 108 rows are retained (fallback for any data
## that doesn't get recoded by the driver patches).
##=============================================================================

suppressPackageStartupMessages({ library(data.table) })

IN_FILE  <- "/users/PUOM0008/crsfaaron/fvs-conus/traits/species_traits.rds"
OUT_FILE <- "/users/PUOM0008/crsfaaron/fvs-conus/traits/species_traits_v2.rds"

traits <- as.data.table(readRDS(IN_FILE))
cat(sprintf("Loaded %d species from %s\n", nrow(traits), IN_FILE))

# ---------------------------------------------------------------------------
# Douglas-fir split
# ---------------------------------------------------------------------------
df_template <- traits[SPCD == 202]
stopifnot(nrow(df_template) == 1)

coastal_df <- copy(df_template)
coastal_df[, SPCD := 2020L]
coastal_df[, SPECIES := "menziesii_coastal"]
coastal_df[, COMMON_NAME := "Coastal Douglas-fir (var. menziesii)"]
coastal_df[, wood_specific_gravity := 0.45]
coastal_df[, shade_tolerance_num := 2.8]
coastal_df[, max_ht_m := 100]
coastal_df[, max_dbh_cm := 200]
coastal_df[, leaf_longevity_months := 60]
coastal_df[, max_ht_fia_p99_m := 100]
coastal_df[, max_dbh_fia_p99_cm := 200]

rocky_df <- copy(df_template)
rocky_df[, SPCD := 2021L]
rocky_df[, SPECIES := "menziesii_glauca"]
rocky_df[, COMMON_NAME := "Rocky Mountain Douglas-fir (var. glauca)"]
rocky_df[, wood_specific_gravity := 0.46]
rocky_df[, shade_tolerance_num := 2.5]
rocky_df[, max_ht_m := 60]
rocky_df[, max_dbh_cm := 130]
rocky_df[, leaf_longevity_months := 72]
rocky_df[, max_ht_fia_p99_m := 60]
rocky_df[, max_dbh_fia_p99_cm := 130]

# ---------------------------------------------------------------------------
# Lodgepole pine split
# ---------------------------------------------------------------------------
lp_template <- traits[SPCD == 108]
stopifnot(nrow(lp_template) == 1)

# Shore pine (var. contorta): coastal, gnarly, much smaller, fully intolerant.
# Wood density slightly higher (slower growth in coastal acid-bog/saltspray sites).
shore_lp <- copy(lp_template)
shore_lp[, SPCD := 1080L]
shore_lp[, SPECIES := "contorta_contorta"]
shore_lp[, COMMON_NAME := "Shore pine (var. contorta)"]
shore_lp[, wood_specific_gravity := 0.40]
shore_lp[, shade_tolerance_num := 1.5]    # very intolerant, slightly more so than interior
shore_lp[, max_ht_m := 25]                # typical mature shore pine
shore_lp[, max_dbh_cm := 50]              # typical mature
shore_lp[, leaf_longevity_months := 60]
shore_lp[, max_ht_fia_p99_m := 25]
shore_lp[, max_dbh_fia_p99_cm := 50]

# Rocky Mountain / Sierra lodgepole (var. latifolia + murrayana): straight bole,
# moderate height, intolerant. Sierra lodgepole is morphologically intermediate
# but operationally lumped here given small Sierra share in the dataset.
rocky_lp <- copy(lp_template)
rocky_lp[, SPCD := 1081L]
rocky_lp[, SPECIES := "contorta_latifolia"]
rocky_lp[, COMMON_NAME := "Rocky Mountain lodgepole pine (var. latifolia)"]
rocky_lp[, wood_specific_gravity := 0.38]
rocky_lp[, shade_tolerance_num := 1.6]
rocky_lp[, max_ht_m := 50]                # typical mature interior lodgepole
rocky_lp[, max_dbh_cm := 80]              # typical mature
rocky_lp[, leaf_longevity_months := 60]
rocky_lp[, max_ht_fia_p99_m := 50]
rocky_lp[, max_dbh_fia_p99_cm := 80]

# ---------------------------------------------------------------------------
# Append all four variety rows
# ---------------------------------------------------------------------------
traits_v2 <- rbind(traits, coastal_df, rocky_df, shore_lp, rocky_lp, fill = TRUE)
cat(sprintf("\nAdded 4 variety rows. New trait file has %d species.\n", nrow(traits_v2)))

cat("\n=== DF variety rows in traits_v2 ===\n")
print(traits_v2[SPCD %in% c(202, 2020, 2021),
                .(SPCD, GENUS, SPECIES, COMMON_NAME,
                  wood_specific_gravity, shade_tolerance_num,
                  max_ht_m, max_dbh_cm, leaf_longevity_months)])

cat("\n=== Lodgepole variety rows in traits_v2 ===\n")
print(traits_v2[SPCD %in% c(108, 1080, 1081),
                .(SPCD, GENUS, SPECIES, COMMON_NAME,
                  wood_specific_gravity, shade_tolerance_num,
                  max_ht_m, max_dbh_cm, leaf_longevity_months)])

saveRDS(traits_v2, OUT_FILE)
cat(sprintf("\nWrote: %s\n", OUT_FILE))
