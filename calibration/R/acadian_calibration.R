# =============================================================================
# calibration/R/acadian_calibration.R
#
# Apply the fvs-modern ACD calibration to the Acadian (AcadianGY / FVS-ACD)
# model through its NATIVE per-tree multiplier inputs.
#
# Why this exists
# ---------------
# FVS-ACD does not honor FVS keyword multipliers (a stress test applying a 5x
# MORTMULT moved 100-year ACD biomass by only ~6%). The Acadian model instead
# reads per-tree dDBH.mult / dHt.mult / mort.mult columns and applies them
# internally. So the ACD calibration must be wired through those columns rather
# than through generate_keywords(). This module maps the ACD
# calibration_multipliers block onto the Acadian tree inputs by species:
#
#   tree SP (Acadian 2-letter code)
#     -> SPCD                (species crosswalk, FIA column)
#     -> FVS species slot    (config categories.species_definitions.FIAJSP)
#     -> dDBH.mult <- dds_multiplier[slot]   (diameter growth)
#        dHt.mult  <- htg_multiplier[slot]   (height increment)
#        mort.mult <- mort_multiplier[slot]  (mortality)
#
# Species without a calibration factor (or not in the crosswalk) default to 1.0.
# =============================================================================

#' Build an Acadian-species lookup of calibration multipliers.
#'
#' @param acd_config_path path to config/calibrated/acd.json
#' @param crosswalk_path path to the Acadian species crosswalk CSV (needs an
#'   Acadian 2-letter code column, e.g. AD_PSP, and an FIA/SPCD column).
#' @return data.frame(SP, dDBH.mult, dHt.mult, mort.mult)
acadian_calibration_lookup <- function(acd_config_path, crosswalk_path) {
  cfg <- jsonlite::fromJSON(acd_config_path, simplifyVector = TRUE)
  cm <- cfg$calibration_multipliers
  if (is.null(cm)) stop("no calibration_multipliers block in ", acd_config_path)
  fia <- suppressWarnings(as.integer(cfg$categories$species_definitions$FIAJSP))

  xw <- utils::read.csv(crosswalk_path, stringsAsFactors = FALSE, check.names = FALSE)
  ad_col <- intersect(c("AD_PSP", "OSM_AD_CmdKey", "AD"), names(xw))
  fia_col <- intersect(c("FIA", "SPCD", "FIAJSP"), names(xw))
  if (!length(ad_col) || !length(fia_col))
    stop("crosswalk missing Acadian-code or FIA/SPCD column: ", crosswalk_path)
  ad <- toupper(trimws(as.character(xw[[ad_col[1]]])))
  spcd <- suppressWarnings(as.integer(xw[[fia_col[1]]]))

  slot_of_spcd <- function(s) {
    w <- which(fia == s)
    if (length(w)) w[1] else NA_integer_
  }
  get <- function(arr, slot) if (!is.na(slot) && slot <= length(arr)) arr[[slot]] else 1.0

  sp_unique <- unique(ad[!is.na(ad) & ad != ""])
  out <- data.frame(SP = sp_unique, dDBH.mult = 1.0, dHt.mult = 1.0,
                    mort.mult = 1.0, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(out))) {
    sc <- spcd[match(out$SP[i], ad)]
    slot <- if (!is.na(sc)) slot_of_spcd(sc) else NA_integer_
    out$dDBH.mult[i] <- get(cm$dds_multiplier, slot)
    out$dHt.mult[i]  <- get(cm$htg_multiplier, slot)
    out$mort.mult[i] <- get(cm$mort_multiplier, slot)
  }
  out
}

#' Set per-tree dDBH.mult / dHt.mult / mort.mult on an Acadian tree data frame.
#'
#' @param tree_df data frame with an SP column (Acadian 2-letter code).
#' @param acd_config_path,crosswalk_path see acadian_calibration_lookup().
#' @param version "calibrated" applies the factors; anything else sets 1.0.
#' @return tree_df with the three multiplier columns populated.
apply_acadian_calibration <- function(tree_df, acd_config_path, crosswalk_path,
                                      version = "calibrated") {
  if (!identical(version, "calibrated")) {
    tree_df$dDBH.mult <- 1.0
    tree_df$dHt.mult <- 1.0
    tree_df$mort.mult <- 1.0
    return(tree_df)
  }
  lut <- acadian_calibration_lookup(acd_config_path, crosswalk_path)
  idx <- match(toupper(trimws(tree_df$SP)), lut$SP)
  tree_df$dDBH.mult <- ifelse(is.na(idx), 1.0, lut$dDBH.mult[idx])
  tree_df$dHt.mult  <- ifelse(is.na(idx), 1.0, lut$dHt.mult[idx])
  tree_df$mort.mult <- ifelse(is.na(idx), 1.0, lut$mort.mult[idx])
  tree_df
}
