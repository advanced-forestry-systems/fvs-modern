## =============================================================================
## 30b_build_htdbh_input.R  (v2)
##
## Alternative to 30_build_conus_dataset.R for the HT-DBH benchmark only.
## Reads standard FIADB tables (no CHANGEdata zips required) and assembles
## cross-sectional tree data with the same column names prepare_htdbh_data()
## expects.
##
## Columns computed here vs CHANGEdata source:
##   BA1      <- sum(pi*(DIA/24)^2 * TPA_UNADJ) per plot (ft^2/ac)
##   BAL1     <- BA of trees with DIA > target tree in same plot (ft^2/ac)
##   SDI1     <- sum(TPA_UNADJ * (DIA/10)^1.605) per plot (Reineke)
##   rd_add   <- SDI1 / SDImax_brms (from ~/SDImax/brms.SDImax.1-27-24.csv)
##   EPA_L3_CODE <- from blup_plot_asymptotes.csv (already has this column)
##   EPA_L1/L2   <- parsed from EPA_L3_CODE (string prefix)
##   cspi     <- weighted composite of BGI_z (0.533) + Asym_z (0.467)
##              (ClimateSI unavailable; weights renormalized)
##
## All primary inputs present on Cardinal:
##   ~/FIA/XX_{TREE,PLOT,COND}.csv   (49 CONUS states, standard FIADB)
##   ~/FIA/asym_agb_analysis/results/stage6_bgi_results.rds
##   ~/FIA/asym_agb_analysis/results/blup_plot_asymptotes.csv  (has EPA_L3_CODE)
##   ~/SDImax/brms.SDImax.1-27-24.csv
##
## Output:
##   calibration/data/conus_remeasurement_pairs.rds
##   calibration/data/conus_htdbh_build_summary.csv
##
## Author: A. Weiskittel, G. Johnson
## Date: 2026-07-09
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glue)
})

## ---- Configuration ----------------------------------------------------------

FIA_DIR    <- path.expand("~/FIA")
SDIMAX_CSV <- path.expand("~/SDImax/brms.SDImax.1-27-24.csv")
OUT_DIR    <- "calibration/data"
BGI_RDS    <- file.path(FIA_DIR, "asym_agb_analysis/results/stage6_bgi_results.rds")
BLUP_CSV   <- file.path(FIA_DIR, "asym_agb_analysis/results/blup_plot_asymptotes.csv")

CONDPROP_MIN <- 0.75
EXCLUDE_AK   <- TRUE
AK_STATECD   <- 2L

## CSPI weights (stage9 analysis: BGI 0.40, Asym 0.35, SI 0.25)
## ClimateSI not available; renormalize to BGI 0.533, Asym 0.467
CSPI_W_BGI  <- 0.40
CSPI_W_ASYM <- 0.35

## TPA_UNADJ approximations for standard DESIGNCD=1 FIA plots
TPA_SUBPLOT  <- 6.018    # trees DIA >= 5" (24 ft radius subplot)
TPA_MICROPLOT <- 74.965  # trees 1" <= DIA < 5" (6.8 ft radius microplot)

## 49 CONUS state abbreviations
STATES <- c("AL","AR","AZ","CA","CO","CT","DE","FL","GA","IA","ID","IL","IN",
            "KS","KY","LA","MA","MD","ME","MI","MN","MO","MS","MT","NC","ND",
            "NE","NH","NJ","NM","NV","NY","OH","OK","OR","PA","RI","SC","SD",
            "TN","TX","UT","VA","VT","WA","WI","WV","WY")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("===============================================================\n")
cat("FVS-CONUS HT-DBH Input Builder v2 (FIADB alternative)\n")
cat("  Start:", format(Sys.time()), "\n")
cat("===============================================================\n\n")

## ---- Step 1: Pre-load support tables --------------------------------------

cat("Step 1: Loading support tables...\n")

## BGI: join via rounded LAT/LON
bgi_raw <- readRDS(BGI_RDS)
bgi_dt  <- as.data.table(bgi_raw$plot_bgi)
bgi_dt[, LAT_r := round(LAT, 2)]
bgi_dt[, LON_r := round(LON, 2)]
bgi_agg <- bgi_dt[, .(bgi = mean(BGI_A_mean, na.rm = TRUE)),
                   by = .(LAT_r, LON_r)]
setkeyv(bgi_agg, c("LAT_r", "LON_r"))
cat("  BGI:     ", format(nrow(bgi_agg), big.mark = ","), "LAT/LON cells\n")
rm(bgi_raw, bgi_dt); gc(verbose = FALSE)

## BLUP: join via plot_key; also read EPA_L3_CODE which is in this CSV
blup_dt <- fread(BLUP_CSV, showProgress = FALSE,
                 select = c("PlotID", "Asym_mean", "EPA_L3_CODE"))
blup_dt[, plot_key := gsub("_", "-", PlotID)]
blup_dt[, PlotID   := NULL]
blup_dt[EPA_L3_CODE == "UNKNOWN", EPA_L3_CODE := NA_character_]
setkeyv(blup_dt, "plot_key")
cat("  BLUP:    ", format(nrow(blup_dt), big.mark = ","), "plots",
    " | EPA_L3 coverage:", sum(!is.na(blup_dt$EPA_L3_CODE)),
    sprintf("(%.0f%%)\n", 100 * mean(!is.na(blup_dt$EPA_L3_CODE))))

## SDImax (brms): join via plot_key
sdi_dt <- fread(SDIMAX_CSV, showProgress = FALSE)
## Detect key columns
key_cols <- c("STATECD","UNITCD","COUNTYCD","PLOT")
sdi_dt[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]
sdi_dt <- sdi_dt[, .(plot_key, SDImax_brms = SDImax.mean)]
setkeyv(sdi_dt, "plot_key")
cat("  SDImax:  ", format(nrow(sdi_dt), big.mark = ","), "plots\n\n")

## ---- Step 2: Read FIADB tables state by state and compute stand metrics ----

cat("Step 2: Reading FIADB + computing BA/BAL/SDI (", length(STATES),
    "states)...\n")

tree_need <- c("CN","PLT_CN","STATECD","UNITCD","COUNTYCD","PLOT",
               "INVYR","CONDID","STATUSCD","SPCD","DIA","HT","HTCD","CR",
               "TPA_UNADJ")
plot_need <- c("CN","STATECD","UNITCD","COUNTYCD","PLOT","LAT","LON")
cond_need <- c("PLT_CN","CONDID","COND_STATUS_CD","CONDPROP_UNADJ")

safe_select <- function(f, want) intersect(want, names(fread(f, nrow = 0L)))

all_trees <- vector("list", length(STATES))

for (i in seq_along(STATES)) {
  st     <- STATES[i]
  tree_f <- file.path(FIA_DIR, glue("{st}_TREE.csv"))
  plot_f <- file.path(FIA_DIR, glue("{st}_PLOT.csv"))
  cond_f <- file.path(FIA_DIR, glue("{st}_COND.csv"))

  if (!all(file.exists(tree_f, plot_f, cond_f))) {
    cat("  ", st, ": missing file(s), skipping\n"); next
  }

  ## -- Read TREE: live trees with measured DIA and HT -----------------------
  sel <- safe_select(tree_f, tree_need)
  tr  <- fread(tree_f, showProgress = FALSE, select = sel)
  tr  <- tr[STATUSCD == 1L & !is.na(DIA) & DIA > 0 & !is.na(HT) & HT > 0]
  if (nrow(tr) == 0L) { cat("  ", st, ": 0 live trees\n"); next }

  ## -- COND filter: forested, sufficient proportion -------------------------
  sel_c <- safe_select(cond_f, cond_need)
  cd    <- fread(cond_f, showProgress = FALSE, select = sel_c)
  cd[, PLT_CN := as.character(PLT_CN)]
  cd_ok <- cd[COND_STATUS_CD == 1L &
                !is.na(CONDPROP_UNADJ) & CONDPROP_UNADJ >= CONDPROP_MIN,
              .(PLT_CN, CONDID, CONDPROP_UNADJ)]
  rm(cd)
  if (nrow(cd_ok) == 0L) { cat("  ", st, ": 0 qualifying conditions\n"); next }

  tr[, PLT_CN := as.character(PLT_CN)]
  setkeyv(cd_ok, c("PLT_CN", "CONDID"))
  setkeyv(tr,    c("PLT_CN", "CONDID"))
  tr <- tr[cd_ok, nomatch = 0L]
  rm(cd_ok)
  if (nrow(tr) == 0L) { cat("  ", st, ": 0 trees after cond filter\n"); next }

  ## -- PLOT: get LAT/LON ----------------------------------------------------
  sel_p <- safe_select(plot_f, plot_need)
  pl    <- fread(plot_f, showProgress = FALSE, select = sel_p)
  pl[, CN := as.character(CN)]
  setnames(pl, "CN", "PLT_CN")
  setkeyv(pl, "PLT_CN")
  tr <- merge(tr, pl[, .(PLT_CN, LAT, LON)], by = "PLT_CN", all.x = FALSE)
  rm(pl)

  ## -- TPA_UNADJ: use column if present, else derive from DIA ---------------
  if ("TPA_UNADJ" %in% names(tr)) {
    tr[is.na(TPA_UNADJ), TPA_UNADJ := fifelse(DIA >= 5.0,
                                               TPA_SUBPLOT, TPA_MICROPLOT)]
  } else {
    tr[, TPA_UNADJ := fifelse(DIA >= 5.0, TPA_SUBPLOT, TPA_MICROPLOT)]
  }

  ## -- Per-tree basal area contribution (ft^2/acre) -------------------------
  tr[, BA_TPA := pi * (DIA / 24)^2 * TPA_UNADJ]

  ## -- Plot-level BA1 and SDI1 ----------------------------------------------
  plot_stats <- tr[, .(
    BA1  = sum(BA_TPA,                    na.rm = TRUE),
    SDI1 = sum(TPA_UNADJ * (DIA/10)^1.605, na.rm = TRUE)
  ), by = PLT_CN]

  tr <- merge(tr, plot_stats, by = "PLT_CN", all.x = TRUE)
  rm(plot_stats)

  ## -- Per-tree BAL1 (BA of trees with strictly larger DIA in same plot) ----
  ## Sort descending DIA within plot; shifted cumsum gives BAL above each tier
  ## Aggregate by (PLT_CN, DIA) first to handle ties correctly.
  dia_agg <- tr[, .(tier_ba = sum(BA_TPA, na.rm = TRUE)),
                by = .(PLT_CN, DIA)]
  setorder(dia_agg, PLT_CN, -DIA)
  dia_agg[, cum_above := shift(cumsum(tier_ba), fill = 0), by = PLT_CN]
  dia_agg[, tier_ba := NULL]
  tr <- merge(tr, dia_agg[, .(PLT_CN, DIA, BAL1 = cum_above)],
              by = c("PLT_CN", "DIA"), all.x = TRUE)
  tr[is.na(BAL1), BAL1 := 0]
  rm(dia_agg)

  ## -- Clean up scratch columns ---------------------------------------------
  tr[, BA_TPA   := NULL]
  tr[, TPA_UNADJ := NULL]

  if (EXCLUDE_AK && "STATECD" %in% names(tr))
    tr <- tr[STATECD != AK_STATECD]

  cat("  ", st, ":", format(nrow(tr), big.mark = ","), "trees |",
      "BA1 med =", round(median(tr$BA1, na.rm = TRUE), 1),
      " SDI1 med =", round(median(tr$SDI1, na.rm = TRUE), 1), "\n")

  all_trees[[i]] <- tr
  rm(tr); gc(verbose = FALSE)
}

cat("\nCombining states...\n")
trees <- rbindlist(all_trees, use.names = TRUE, fill = TRUE)
rm(all_trees); gc(verbose = FALSE)
cat("  Total:", format(nrow(trees), big.mark = ","), "trees,",
    uniqueN(trees$STATECD), "states,",
    format(uniqueN(trees$PLT_CN), big.mark = ","), "plots\n\n")

## ---- Step 3: Plot-level joins (BGI, BLUP+EPA, SDImax) ---------------------

cat("Step 3: Joining plot-level attributes...\n")

## 3a. BGI by rounded LAT/LON
trees[, LAT_r := round(LAT, 2)]
trees[, LON_r := round(LON, 2)]
trees <- merge(trees, bgi_agg, by = c("LAT_r", "LON_r"), all.x = TRUE)
trees[, c("LAT_r", "LON_r") := NULL]
cat("  BGI coverage:     ", sprintf("%.1f%%\n",
    100 * mean(!is.na(trees$bgi))))
rm(bgi_agg); gc(verbose = FALSE)

## 3b. BLUP asymptotes + EPA_L3_CODE by plot_key
trees[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]
trees <- merge(trees, blup_dt, by = "plot_key", all.x = TRUE)
setnames(trees, "Asym_mean", "max_biomass")
cat("  BLUP coverage:    ", sprintf("%.1f%%\n",
    100 * mean(!is.na(trees$max_biomass))))
cat("  EPA_L3 coverage:  ", sprintf("%.1f%%\n",
    100 * mean(!is.na(trees$EPA_L3_CODE))))
rm(blup_dt); gc(verbose = FALSE)

## 3c. SDImax by plot_key -> rd_add
trees <- merge(trees, sdi_dt, by = "plot_key", all.x = TRUE)
trees[, rd_add := fifelse(
  !is.na(SDImax_brms) & SDImax_brms > 0 & !is.na(SDI1) & SDI1 >= 0,
  SDI1 / SDImax_brms,
  NA_real_
)]
cat("  SDImax coverage:  ", sprintf("%.1f%%\n",
    100 * mean(!is.na(trees$SDImax_brms))))
cat("  rd_add coverage:  ", sprintf("%.1f%%", 100 * mean(!is.na(trees$rd_add))))
cat("  | median rd_add:", round(median(trees$rd_add, na.rm=TRUE), 3), "\n")
rm(sdi_dt); gc(verbose = FALSE)

## ---- Step 4: EPA L1/L2 from L3 code, unit conversion, CSPI ---------------

cat("\nStep 4: EPA hierarchy, unit conversion, CSPI...\n")

## EPA L3 code format: "5.2.3" -> L1 = "5", L2 = "5.2"
trees[, EPA_L1_CODE := sub("^([^.]+).*$",         "\\1", EPA_L3_CODE)]
trees[, EPA_L2_CODE := sub("^([^.]+\\.[^.]+).*$", "\\1", EPA_L3_CODE)]
cat("  EPA L1 unique:", uniqueN(trees$EPA_L1_CODE[!is.na(trees$EPA_L1_CODE)]),
    "| L2:", uniqueN(trees$EPA_L2_CODE[!is.na(trees$EPA_L2_CODE)]),
    "| L3:", uniqueN(trees$EPA_L3_CODE[!is.na(trees$EPA_L3_CODE)]), "\n")

## Unit conversion: DIA inches -> cm, HT feet -> m
trees[, DBH1 := DIA * 2.54]
trees[, HT1  := HT  * 0.3048]

## Rename to prepare_htdbh_data() expected names
setnames(trees, c("STATUSCD","HTCD"), c("STATUS1","HTCD1"), skip_absent = TRUE)
if ("CR" %in% names(trees)) {
  setnames(trees, "CR", "CR1")
  trees[, CR1 := CR1 / 100]
}

## CSPI: z-score composite (BGI 0.533, Asym 0.467; ClimateSI unavailable)
w_bgi  <- CSPI_W_BGI  / (CSPI_W_BGI + CSPI_W_ASYM)
w_asym <- CSPI_W_ASYM / (CSPI_W_BGI + CSPI_W_ASYM)

if (any(!is.na(trees$bgi)))
  trees[, bgi_z := (bgi - mean(bgi, na.rm=TRUE)) / sd(bgi, na.rm=TRUE)]  else
  trees[, bgi_z := NA_real_]

if (any(!is.na(trees$max_biomass)))
  trees[, asym_z := (max_biomass - mean(max_biomass, na.rm=TRUE)) /
                    sd(max_biomass, na.rm=TRUE)]  else
  trees[, asym_z := NA_real_]

trees[, cspi := {
  hb <- !is.na(bgi_z); ha <- !is.na(asym_z)
  fifelse(hb & ha, w_bgi * bgi_z + w_asym * asym_z,
  fifelse(hb, bgi_z, fifelse(ha, asym_z, NA_real_)))
}]
trees[, c("bgi_z","asym_z") := NULL]

cat("  CSPI coverage:   ", sprintf("%.1f%%\n", 100*mean(!is.na(trees$cspi))))
cat("  DBH1 range: [",  round(min(trees$DBH1, na.rm=TRUE),1),
    ",", round(max(trees$DBH1, na.rm=TRUE),1), "] cm\n")
cat("  HT1 range:  [",  round(min(trees$HT1,  na.rm=TRUE),2),
    ",", round(max(trees$HT1,  na.rm=TRUE),2), "] m\n")
cat("  BA1 range:  [",  round(min(trees$BA1,  na.rm=TRUE),1),
    ",", round(max(trees$BA1,  na.rm=TRUE),1), "] ft^2/ac\n")
cat("  BAL1 range: [",  round(min(trees$BAL1, na.rm=TRUE),1),
    ",", round(max(trees$BAL1, na.rm=TRUE),1), "] ft^2/ac\n\n")

## ---- Step 5: Save output --------------------------------------------------

cat("Step 5: Saving output...\n")

## Drop raw FIA columns no longer needed
drop_cols <- c("DIA","HT","CN","SDImax_brms","SDI1")
trees[, (intersect(drop_cols, names(trees))) := NULL]

## Per-state summary
summ <- trees[, .(
  n_trees    = .N,
  n_plots    = uniqueN(PLT_CN),
  n_spcd     = uniqueN(SPCD),
  pct_htcd1  = if ("HTCD1" %in% names(.SD))
                 round(100 * mean(HTCD1 == 1L, na.rm=TRUE), 1) else NA_real_,
  pct_cspi   = round(100 * mean(!is.na(cspi)), 1),
  pct_ba1    = round(100 * mean(!is.na(BA1) & BA1 > 0), 1),
  pct_rdadd  = round(100 * mean(!is.na(rd_add)), 1),
  pct_epa    = round(100 * mean(!is.na(EPA_L3_CODE)), 1),
  dbh1_med   = round(median(DBH1, na.rm=TRUE), 2),
  ht1_med    = round(median(HT1,  na.rm=TRUE), 2),
  ba1_med    = round(median(BA1,  na.rm=TRUE), 1),
  rd_med     = round(median(rd_add, na.rm=TRUE), 3)
), by = STATECD]

fwrite(summ, file.path(OUT_DIR, "conus_htdbh_build_summary.csv"))
cat("  Summary ->", file.path(OUT_DIR, "conus_htdbh_build_summary.csv"), "\n")

saveRDS(trees, file.path(OUT_DIR, "conus_remeasurement_pairs.rds"))
cat("  Pairs   ->", file.path(OUT_DIR, "conus_remeasurement_pairs.rds"), "\n")
cat("  Rows:", format(nrow(trees), big.mark=","),
    "| Cols:", ncol(trees), "\n")
cat("  Columns:", paste(sort(names(trees)), collapse=", "), "\n\n")
cat("=== BUILD COMPLETE:", format(Sys.time()), "===\n")
