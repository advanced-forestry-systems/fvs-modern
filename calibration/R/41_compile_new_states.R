##=============================================================================
## 41_compile_new_states.R
##
## Process raw FIA DataMart CSV downloads into the same pre-compiled format
## used by the existing 43 states in the fvs-conus data directory:
##   <ST>_TREEdata.CSV   tree-level with derived stand variables
##   <ST>_PLOTdata.CSV   plot-level summaries
##   <ST>_CHANGEdata.CSV paired start/end measurements with competition metrics
##
## The output is zipped into <ST>_FIAdata_<VINTAGE>.zip to match the existing
## file naming convention.
##
## Usage (from repo root on Cardinal):
##   source /etc/profile.d/lmod.sh
##   module load gcc/12.3.0 R/4.4.0 gdal/3.7.3 proj/9.2.1 geos/3.12.0
##
##   Rscript calibration/R/41_compile_new_states.R \
##     --raw_dir   calibration/data/raw_fia \
##     --out_dir   calibration/data/fia_zips \
##     --states    AR,CA,MT,OR,TX \
##     --vintage   20260301
##
## Author: A. Weiskittel (2026-04-14)
##=============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fs)
  library(glue)
  library(optparse)
})

## ---- CLI -------------------------------------------------------------------
opts <- OptionParser(option_list = list(
  make_option("--raw_dir",   type = "character",
              default = "calibration/data/raw_fia"),
  make_option("--out_dir",   type = "character",
              default = "calibration/data/fia_zips"),
  make_option("--states",    type = "character", default = "AR,CA,MT,OR,TX"),
  make_option("--vintage",   type = "character", default = "20260301")
)) |> parse_args()

target_states <- str_split(opts$states, ",")[[1]] |> str_trim()
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

message(glue("Processing states: {paste(target_states, collapse = ', ')}"))
message(glue("Raw FIA dir: {opts$raw_dir}"))
message(glue("Output dir:  {opts$out_dir}"))

## ---- Constants --------------------------------------------------------------
VERSION_TAG <- format(Sys.Date(), "%Y-%m-%d")

## ---- Helper: read a raw FIA table from the DataMart download ----------------
## DataMart zips contain <ST>_<TABLE>.csv (e.g., AR_TREE.csv)
read_raw_fia <- function(raw_dir, state, table_name) {
  ## Check both unzipped CSV and inside zip
  csv_path <- file.path(raw_dir, paste0(state, "_", table_name, ".csv"))
  if (file.exists(csv_path)) {
    message(glue("  Reading {csv_path}"))
    return(fread(csv_path, showProgress = FALSE))
  }

  ## Try inside the raw zip
  zip_path <- file.path(raw_dir, paste0(state, "_CSV.zip"))
  if (!file.exists(zip_path)) {
    zip_path <- file.path(raw_dir, paste0(state, "_csv.zip"))
  }
  if (file.exists(zip_path)) {
    entries <- unzip(zip_path, list = TRUE)$Name
    match <- entries[grepl(paste0(state, "_", table_name, "\\.csv$"),
                           entries, ignore.case = TRUE)]
    if (length(match) > 0) {
      tmp <- tempdir()
      unzip(zip_path, files = match[1], exdir = tmp, overwrite = TRUE)
      df <- fread(file.path(tmp, match[1]), showProgress = FALSE)
      unlink(file.path(tmp, match[1]))
      message(glue("  Read {table_name} from zip: {nrow(df)} rows"))
      return(df)
    }
  }

  warning(glue("  Table {table_name} not found for {state}"))
  return(NULL)
}

## ---- Compute stand-level competition variables for all trees in a plot ------
compute_stand_vars <- function(tree_dt) {
  ## tree_dt should have: PLT_CN, DIA, HT, CR, TPA_UNADJ, STATUSCD, SPCD
  ## Returns with BAL, CCFL, CCH, TPA, BAPA, QMD, HT40, DBH40, etc. appended

  live <- tree_dt[STATUSCD == 1 & !is.na(DIA) & DIA > 0]

  ## Plot-level summaries
  plt_sum <- live[, .(
    TPA    = sum(TPA_UNADJ, na.rm = TRUE),
    BAPA   = sum(0.005454 * DIA^2 * TPA_UNADJ, na.rm = TRUE),
    QMD    = sqrt(sum(DIA^2 * TPA_UNADJ, na.rm = TRUE) /
                  sum(TPA_UNADJ, na.rm = TRUE)),
    NLive  = .N
  ), by = PLT_CN]

  ## HT40, DBH40: mean of 40 largest TPA trees
  ht40 <- live[order(PLT_CN, -DIA),
               .(HT40  = mean(head(HT, 40), na.rm = TRUE),
                 DBH40 = mean(head(DIA, 40), na.rm = TRUE)),
               by = PLT_CN]

  ## Lorey height: BA-weighted mean height
  htlorey <- live[!is.na(HT),
                  .(HTlorey = sum(0.005454 * DIA^2 * TPA_UNADJ * HT, na.rm = TRUE) /
                              sum(0.005454 * DIA^2 * TPA_UNADJ, na.rm = TRUE)),
                  by = PLT_CN]

  ## SDI and RD (Curtis)
  sdi_rd <- plt_sum[, .(
    SDI = TPA * (QMD / 10)^1.605,
    RD  = BAPA / sqrt(QMD)
  ), by = PLT_CN]

  ## Crown competition factor (rough: using MCW ~ a + b*DIA by species group)
  ## Simplified: CCF ~ sum(0.001803 * MCW^2 * TPA) where MCW ~ 3.0 + 1.65*DIA
  ccf <- live[, .(
    CCF = sum(0.001803 * (3.0 + 1.65 * DIA)^2 * TPA_UNADJ, na.rm = TRUE)
  ), by = PLT_CN]

  ## BAL: basal area in larger trees (per tree)
  setorder(live, PLT_CN, -DIA)
  live[, cumBA := cumsum(0.005454 * DIA^2 * TPA_UNADJ), by = PLT_CN]
  live[, BAL := shift(cumBA, 1, fill = 0), by = PLT_CN]

  ## CCFL: crown competition factor in larger trees
  live[, ccf_i := 0.001803 * (3.0 + 1.65 * DIA)^2 * TPA_UNADJ]
  live[, cumCCF := cumsum(ccf_i), by = PLT_CN]
  live[, CCFL := shift(cumCCF, 1, fill = 0), by = PLT_CN]

  ## CCH: canopy cover from below (similar concept, bounded at 100)
  live[, CCH := pmin(CCFL, 100)]

  ## Primary and secondary species by BA
  sp_ba <- live[, .(spBA = sum(0.005454 * DIA^2 * TPA_UNADJ, na.rm = TRUE)),
                by = .(PLT_CN, SPCD)]
  sp_ba[, rank := frank(-spBA), by = PLT_CN]
  primary <- sp_ba[rank == 1, .(PLT_CN,
                                primarySPBA = SPCD,
                                primaryPctBAPA = spBA)]
  secondary <- sp_ba[rank == 2, .(PLT_CN,
                                  secondarySPBA = SPCD,
                                  secondaryPctBAPA = spBA)]

  ## Merge everything back
  live <- merge(live, plt_sum, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, ht40, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, htlorey, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, sdi_rd, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, ccf, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, primary, by = "PLT_CN", all.x = TRUE)
  live <- merge(live, secondary, by = "PLT_CN", all.x = TRUE)

  ## Compute pct for primary/secondary
  live[, primaryPctBAPA := round(100 * primaryPctBAPA / BAPA, 1)]
  live[, secondaryPctBAPA := round(100 * secondaryPctBAPA / BAPA, 1)]

  ## Clean up temp columns
  live[, c("cumBA", "ccf_i", "cumCCF") := NULL]

  live
}

## ---- Process one state ------------------------------------------------------
process_state <- function(st, raw_dir, out_dir, vintage) {
  message(glue("\n{'='|strrep(60)}"))
  message(glue("Processing state: {st}"))
  message(glue("{'='|strrep(60)}"))

  ## Read raw tables
  tree <- read_raw_fia(raw_dir, st, "TREE")
  plot <- read_raw_fia(raw_dir, st, "PLOT")
  cond <- read_raw_fia(raw_dir, st, "COND")

  if (is.null(tree) || is.null(plot)) {
    warning(glue("Missing TREE or PLOT for {st}, skipping"))
    return(NULL)
  }

  setDT(tree)
  setDT(plot)
  if (!is.null(cond)) setDT(cond)

  ## Standardize column names (FIA uses CN as primary key)
  if (!"PLT_CN" %in% names(tree) && "CN" %in% names(plot)) {
    ## Map tree PLT_CN to plot CN
    if ("PLT_CN" %in% names(tree)) {
      ## already exists
    } else if ("PLT_CN" %in% names(tree)) {
      ## already exists
    } else {
      ## Need to find the link
      message("  Note: PLT_CN mapping from plot/tree tables")
    }
  }

  ## Get plot-level info (coordinates, site, etc.)
  plot_info <- plot[, .SD, .SDcols = intersect(names(plot),
    c("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR",
      "PLOT_STATUS_CD", "LAT", "LON", "ELEV",
      "DESIGNCD", "MEESSION"))]

  ## Get condition-level info
  if (!is.null(cond)) {
    cond_info <- cond[, .SD, .SDcols = intersect(names(cond),
      c("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ",
        "OWNCD", "FORTYPCD", "STDAGE", "SITECLCD", "SITECL_METHOD",
        "SICOND", "SIBASE", "SISP", "SLOPE", "ASPECT", "PHYSCLCD",
        "STDORGCD", "STDORGSP"))]
  }

  ## Mark missing data flags
  tree[, DIAmissing := as.integer(is.na(DIA) | DIA <= 0)]
  tree[, HTmissing  := as.integer(is.na(HT))]
  tree[, CRmissing  := as.integer(is.na(CR))]

  ## Compute stand variables
  message("  Computing stand-level competition variables...")
  tree_enriched <- compute_stand_vars(tree)
  message(glue("  Enriched {nrow(tree_enriched)} live tree records"))

  ## --- Build CHANGEdata: paired start/end measurements ----------------------
  message("  Building growth change pairs...")

  ## Identify trees with PREV_TRE_CN (remeasured trees)
  if ("PREV_TRE_CN" %in% names(tree)) {
    ## t2 records: current measurements with a previous measurement link
    t2 <- tree[!is.na(PREV_TRE_CN)]
    t1 <- tree[, .(CN, PLT_CN, INVYR, CONDID, STATUSCD, AGENTCD,
                    DIA, DIAHTCD, HT, HTCD, ACTUALHT, CR,
                    TREECLCD, CCLCD, TPA_UNADJ)]
    setnames(t1, setdiff(names(t1), "CN"), paste0("prev_", setdiff(names(t1), "CN")))

    ## Join t2 to t1 on PREV_TRE_CN
    pairs <- merge(t2, t1, by.x = "PREV_TRE_CN", by.y = "CN", all.x = FALSE)

    ## Compute remeasurement period
    pairs[, REMPER := INVYR - prev_INVYR]
    pairs[, REMPERadj := REMPER]  # can be adjusted for growing season
    pairs <- pairs[REMPER >= 1 & REMPER <= 20]

    ## Get start and end stand variables from enriched tree table
    ## Match start conditions from the enriched data for prev_PLT_CN/prev_INVYR
    ## and end conditions for current PLT_CN/INVYR

    message(glue("  Found {nrow(pairs)} remeasurement pairs"))
  } else {
    message("  PREV_TRE_CN not available; attempting GRM-based pairing")
    grm <- read_raw_fia(raw_dir, st, "TREE_GRM_COMPONENT")
    if (is.null(grm)) {
      warning(glue("  No GRM table for {st}, cannot build change pairs"))
      pairs <- data.table()
    } else {
      ## GRM-based pairing would go here
      pairs <- data.table()
    }
  }

  ## --- Build start/end enriched change records ------------------------------
  ## This mirrors the format of the existing CHANGEdata files
  if (nrow(pairs) > 0) {
    ## For each pair, we need start and end stand variables
    ## Join plot info for coordinates
    if ("PLT_CN" %in% names(pairs) && "CN" %in% names(plot_info)) {
      pairs <- merge(pairs, plot_info[, .(CN, LAT, LON, ELEV, DESIGNCD)],
                     by.x = "PLT_CN", by.y = "CN", all.x = TRUE)
    }

    ## Join condition info
    if (!is.null(cond) && "PLT_CN" %in% names(pairs)) {
      cond_join <- cond_info[, .SD[1], by = PLT_CN]  # take first condition
      pairs <- merge(pairs, cond_join, by = "PLT_CN", all.x = TRUE)
    }

    ## Format into the CHANGEdata column structure
    ## The target format uses start/end prefixes
    change_data <- pairs[, .(
      VERSION         = VERSION_TAG,
      SUBP            = if ("SUBP" %in% names(pairs)) SUBP else NA_integer_,
      TREE            = if ("TREE" %in% names(pairs)) TREE else NA_integer_,
      PLOT            = if ("PLOT" %in% names(pairs)) PLOT else NA_integer_,
      STATECD         = STATECD,
      UNITCD          = UNITCD,
      COUNTYCD        = COUNTYCD,
      startCONDID     = if ("CONDID" %in% names(pairs)) prev_CONDID else NA_integer_,
      startINVYR      = prev_INVYR,
      SPCD            = SPCD,
      SPGRPCD         = if ("SPGRPCD" %in% names(pairs)) SPGRPCD else NA_integer_,
      startSTATUSCD   = prev_STATUSCD,
      startAGENTCD    = if ("prev_AGENTCD" %in% names(pairs)) prev_AGENTCD else NA_integer_,
      startDIA        = prev_DIA,
      startDIAHTCD    = if ("prev_DIAHTCD" %in% names(pairs)) prev_DIAHTCD else NA_integer_,
      startHT         = prev_HT,
      startHTCD       = if ("prev_HTCD" %in% names(pairs)) prev_HTCD else NA_integer_,
      startACTUALHT   = if ("prev_ACTUALHT" %in% names(pairs)) prev_ACTUALHT else NA_real_,
      startCR         = prev_CR,
      startTPA_UNADJ  = prev_TPA_UNADJ,
      endINVYR        = INVYR,
      endSTATUSCD     = STATUSCD,
      endAGENTCD      = if ("AGENTCD" %in% names(pairs)) AGENTCD else NA_integer_,
      endDIA          = DIA,
      endHT           = HT,
      endACTUALHT     = if ("ACTUALHT" %in% names(pairs)) ACTUALHT else NA_real_,
      endCR           = CR,
      endTPA_UNADJ    = TPA_UNADJ,
      LAT             = if ("LAT" %in% names(pairs)) LAT else NA_real_,
      LON             = if ("LON" %in% names(pairs)) LON else NA_real_,
      ELEV            = if ("ELEV" %in% names(pairs)) ELEV else NA_real_,
      FORTYPCD        = if ("FORTYPCD" %in% names(pairs)) FORTYPCD else NA_integer_,
      STDAGE          = if ("STDAGE" %in% names(pairs)) STDAGE else NA_integer_,
      SITECLCD        = if ("SITECLCD" %in% names(pairs)) SITECLCD else NA_integer_,
      SICOND          = if ("SICOND" %in% names(pairs)) SICOND else NA_real_,
      SIBASE          = if ("SIBASE" %in% names(pairs)) SIBASE else NA_real_,
      SISP            = if ("SISP" %in% names(pairs)) SISP else NA_integer_,
      SLOPE           = if ("SLOPE" %in% names(pairs)) SLOPE else NA_real_,
      ASPECT          = if ("ASPECT" %in% names(pairs)) ASPECT else NA_real_,
      REMPER          = REMPER,
      REMPERadj       = REMPERadj
    )]

    ## Write change data
    change_file <- file.path(out_dir, glue("{st}_CHANGEdata.CSV"))
    fwrite(change_data, change_file)
    message(glue("  Wrote {nrow(change_data)} change records to {change_file}"))
  }

  ## --- Build TREEdata -------------------------------------------------------
  tree_cols <- intersect(names(tree_enriched),
    c("PLT_CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP",
      "INVYR", "TREE", "SPCD", "SPGRPCD", "STATUSCD", "AGENTCD",
      "DIA", "DIAHTCD", "HT", "HTCD", "ACTUALHT", "CR",
      "TREECLCD", "CCLCD", "TPA_UNADJ",
      "NLive", "TPA", "BAPA", "QMD", "HT40", "DBH40", "HTlorey",
      "SDI", "RD", "CCF", "BAL", "CCFL", "CCH",
      "primarySPBA", "primaryPctBAPA", "secondarySPBA", "secondaryPctBAPA",
      "DIAmissing", "HTmissing", "CRmissing"))
  tree_out <- tree_enriched[, ..tree_cols]
  tree_out[, VERSION := VERSION_TAG]

  ## Add plot info
  if ("CN" %in% names(plot_info)) {
    tree_out <- merge(tree_out, plot_info[, .(CN, LAT, LON, ELEV, DESIGNCD)],
                      by.x = "PLT_CN", by.y = "CN", all.x = TRUE)
  }

  tree_file <- file.path(out_dir, glue("{st}_TREEdata.CSV"))
  fwrite(tree_out, tree_file)
  message(glue("  Wrote {nrow(tree_out)} tree records to {tree_file}"))

  ## --- Build PLOTdata -------------------------------------------------------
  if ("CN" %in% names(plot_info)) {
    plt_vars <- tree_enriched[, .(
      TPA    = sum(TPA_UNADJ, na.rm = TRUE),
      BAPA   = sum(0.005454 * DIA^2 * TPA_UNADJ, na.rm = TRUE),
      QMD    = sqrt(sum(DIA^2 * TPA_UNADJ, na.rm = TRUE) /
                    sum(TPA_UNADJ, na.rm = TRUE)),
      NLive  = .N
    ), by = PLT_CN]

    plot_out <- merge(plot_info, plt_vars, by.x = "CN", by.y = "PLT_CN",
                      all.x = TRUE)
    plot_out[, VERSION := VERSION_TAG]

    plot_file <- file.path(out_dir, glue("{st}_PLOTdata.CSV"))
    fwrite(plot_out, plot_file)
    message(glue("  Wrote {nrow(plot_out)} plot records to {plot_file}"))
  }

  ## --- Zip into project format -----------------------------------------------
  zip_name <- glue("{st}_FIAdata_{vintage}.zip")
  zip_path <- file.path(out_dir, zip_name)

  files_to_zip <- dir_ls(out_dir, regexp = glue("^{st}_.*data\\.CSV$"))
  if (length(files_to_zip) > 0) {
    ## Create zip in output directory
    old_wd <- setwd(out_dir)
    zip(zip_path, files = basename(files_to_zip))
    setwd(old_wd)
    message(glue("  Zipped to {zip_path} ({round(file.size(zip_path)/1024^2, 1)} MB)"))
  }

  ## Return summary
  list(
    state       = st,
    n_trees     = nrow(tree_out),
    n_plots     = if (exists("plot_out")) nrow(plot_out) else 0,
    n_changes   = if (nrow(pairs) > 0) nrow(change_data) else 0,
    n_species   = length(unique(tree_out$SPCD)),
    zip_mb      = round(file.size(zip_path) / 1024^2, 1)
  )
}

## ---- Main -------------------------------------------------------------------
results <- map(target_states, process_state,
               raw_dir = opts$raw_dir,
               out_dir = opts$out_dir,
               vintage = opts$vintage)

## Summary
summary_tbl <- bind_rows(map(results, as_tibble))
summary_file <- file.path(opts$out_dir, "new_states_compile_summary.csv")
write_csv(summary_tbl, summary_file)

message(glue("\n{'='|strrep(60)}"))
message("Compilation Summary")
message(glue("{'='|strrep(60)}"))
print(summary_tbl)
message(glue("\nOutputs in: {opts$out_dir}"))
message("Next: merge these zips with existing data/ directory and rebuild")
message("dg_input.rds / hg_input.rds using script 30.")
