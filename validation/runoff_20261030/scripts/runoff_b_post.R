## runoff_b_post.R OUTROOT VAR ARM MODE BIN
## Adds dominant (top) height to the projector metrics from its own tree lists, by the FVS ATTOPHT rule
## (TPA-weighted mean HT of the 40 largest-DBH trees per acre), so B's HT and form-factor volume are
## defined exactly as for the engines and the observations. long: writes a projector-schema metrics file
## with HT_DOM_M (read by bakuzis_stress_v2.R --source=projector). near: writes the shared near schema.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); ROOT <- a[1]; VAR <- toupper(a[2]); ARM <- a[3]; MODE <- a[4]; BIN <- a[5]; v <- tolower(VAR)
VER <- if (length(a) >= 6) a[6] else "v2"; ROOT <- file.path(ROOT, VER)   # projector version (v2 or v6)
lab <- function(cfg) sub("^conus_", paste0("conus_", VER, "_"), cfg)
topht <- function(d, h, w) { ok <- is.finite(d) & is.finite(h) & is.finite(w) & w > 0; d <- d[ok]; h <- h[ok]; w <- w[ok]
  if (!length(d)) return(NA_real_); o <- order(d, decreasing = TRUE); d <- d[o]; h <- h[o]; w <- w[o]
  cw <- cumsum(w); take <- pmin(w, pmax(40 - (cw - w), 0)); if (sum(take) <= 0) return(NA_real_); sum(h * take) / sum(take) }
dirs <- if (MODE == "long") file.path(ROOT, "long", ARM, v) else grep("/k[0-9]+$", list.dirs(file.path(ROOT, "near", ARM, v), recursive = FALSE), value = TRUE)
res <- rbindlist(lapply(dirs, function(dd) {
  mf <- list.files(dd, "_metrics\\.csv$", full.names = TRUE); tf <- list.files(dd, "_treelists\\.csv$", full.names = TRUE)
  if (!length(mf) || !length(tf)) return(NULL)
  m <- fread(mf[1], colClasses = list(character = "STAND_CN"))
  tl <- fread(tf[1], colClasses = list(character = "STAND_CN"), select = c("STAND_CN", "PROJ_YEAR", "DBH_IN", "HT_M", "TPA"))
  h <- tl[, .(HT_DOM_M = topht(DBH_IN, HT_M, TPA)), by = .(STAND_CN, PROJ_YEAR)]
  m <- merge(m, h, by = c("STAND_CN", "PROJ_YEAR"), all.x = TRUE); m[, CONFIG := lab(CONFIG)]; m }), fill = TRUE)
si <- fread(file.path(BIN, paste0("standinit_", VAR, ".csv")), colClasses = list(character = c("STAND_CN", "STAND_ID")))
cat(sprintf("%s %s %s: standinit stands %d, projected stands %d\n", VAR, ARM, MODE, nrow(si), uniqueN(res$STAND_CN)))
if (MODE == "long") {
  ## Same rule as the engine writer (runoff_engine_runs.py / conus_catalog_runs.py drop TPH <= 0 rows): a stand whose
  ## last tree has died ends its trajectory. With ingrowth off this happens in B; the count is reported.
  dead <- res[!is.finite(BA_FT2AC) | !is.finite(TPH) | TPH <= 0]
  cat(sprintf("rows with no live trees dropped: %d; stands with no live trees by year 100: %d of %d\n", nrow(dead),
      uniqueN(dead[PROJ_YEAR <= 100]$STAND_CN), uniqueN(res$STAND_CN)))
  res <- res[is.finite(BA_FT2AC) & is.finite(TPH) & TPH > 0]
  fwrite(res, file.path(ROOT, "long", sprintf("runoff_long_b_%s_%s_%s_metrics.csv", VER, v, ARM)))
} else {
  b <- res[PROJ_YEAR == 0, .(STAND_CN, ba0_ft2ac = BA_FT2AC, BA0 = BA_FT2AC * 0.2295684, TPH0 = TPH, HT0 = HT_DOM_M)]
  e <- res[PROJ_YEAR > 0, .(STAND_CN, years_run = PROJ_YEAR, config = CONFIG, BA = BA_FT2AC * 0.2295684, TPH = TPH, HT = HT_DOM_M)]
  o <- merge(merge(b, e, by = "STAND_CN"), si[, .(STAND_CN, plot_key, k_cycle = K_CYCLE, balive = BALIVE)], by = "STAND_CN", all.x = TRUE)
  o[, `:=`(variant = VAR, traj = STAND_CN, arm = config,
           seed_ok = is.finite(balive) & balive > 0 & (abs(ba0_ft2ac - balive) / balive <= 0.20 | abs(ba0_ft2ac - balive) <= 1.0),
           QMD0 = sqrt(BA0 / TPH0 / (pi / 40000)), QMD = sqrt(BA / TPH / (pi / 40000)), VOLC = 0.40 * BA * HT)]
  o[, STAND_CN := NULL]
  fwrite(o, file.path(ROOT, "near", sprintf("runoff_near_%s_%s_%s.csv", VER, v, ARM)))
}
