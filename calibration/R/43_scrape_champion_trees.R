#!/usr/bin/env Rscript
# Scrape NCTP 2024 Register and match to FIA SPCD
# Usage: Rscript scripts/43_scrape_champion_trees.R

library(rvest)
library(dplyr)
library(readr)

# Download CSV from NCTP
url <- "http://nationalchampiontree.org/wp-content/uploads/sites/270/2025/02/NCTP-2024-Register-downloadable-measurement-data.csv"
nctp <- read_csv(url, show_col_types = FALSE)

# Load FIA REF_SPECIES (user must provide this file locally)
ref_species <- read_csv("calibration/traits/REF_SPECIES.csv",
                        col_types = "iccc")

# Extract genus + species from NCTP scientific name, convert measurements
champion_data <- nctp %>%
  separate(scientific_name, c("genus", "species"), sep = " ") %>%
  mutate(
    height_m = height_ft * 0.3048,
    dbh_cm = (circumference_in / pi) * 2.54,
    crown_width_m = crown_spread_ft * 0.3048
  ) %>%
  inner_join(ref_species, by = c("genus" = "GENUS", "species" = "SPECIES")) %>%
  select(SPCD, common_name, genus, species, height_m, dbh_cm, crown_width_m) %>%
  mutate(SOURCE = "NCTP 2024 Register", YEAR = 2024) %>%
  arrange(SPCD)

# Write output
write_csv(champion_data, "calibration/traits/champion_trees.csv")
cat(sprintf("Updated champion_trees.csv with %d species\n", nrow(champion_data)))
