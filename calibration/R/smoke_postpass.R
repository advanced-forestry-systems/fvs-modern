#!/usr/bin/env Rscript
# Smoke test: per-stratum ACD post-pass in 19_fia_benchmark_engine.R
# Builds a minimal synthetic validation_data and runs only the post-pass block
# to catch bugs before the 2-hour A/B chain hits them.

suppressMessages(library(data.table))

cat("Smoke test: per-stratum ACD post-pass\n")
cat(strrep("-", 60), "\n", sep = "")

project_root <- "/users/PUOM0008/crsfaaron/fvs-modern-acdbridge"
factors_path <- file.path(project_root,
  "calibration/analysis/acd_stand_level_2026-05-16/acd_calibration_factors.csv")

if (!file.exists(factors_path)) {
  stop("factors CSV missing: ", factors_path)
}

# Build a synthetic validation_data with edge cases:
#  - 10 ACD rows spanning all 4 BA classes, 3 SI terciles, 4 intervals, several FT groups
#  - 5 NE rows (must NOT be touched by post-pass)
#  - 1 row with NA FORTYPCD (should map to "nonstocked-other")
#  - 1 row with NA SICOND (should still get classified)
set.seed(123)
n_acd <- 10
n_ne  <- 5
n_tot <- n_acd + n_ne

validation_data <- data.table(
  VARIANT       = c(rep("ACD", n_acd), rep("NE", n_ne)),
  BA_t1         = c(30, 60, 120, 180, 80, 110, 40, 160, 95, 70,    100, 120, 140, 60, 90),
  SICOND        = c(50, 60, 75, 90, NA, 65, 55, 85, 70, 60,         60, 65, 70, 55, 65),
  FORTYPCD      = c(121L, 801L, 901L, 401L, NA_integer_, 501L, 105L, 121L, 803L, 905L,
                    121L, 801L, 401L, 105L, 121L),
  interval_years = c(5L, 5L, 6L, 7L, 4L, 6L, 5L, 7L, 5L, 5L,        5L, 6L, 5L, 5L, 5L),
  BA_pred_calib  = rep(100.0, n_tot),
  TPA_pred_calib = rep(300.0, n_tot),
  QMD_pred_calib = rep(7.0,   n_tot),
  HT_top_calib   = rep(60.0,  n_tot),
  BA_pred_default = rep(100.0, n_tot)   # filler so "Remove NAs" line works
)
validation_data_save <- copy(validation_data)

# Apply the per-stratum block from 19_fia_benchmark_engine.R
# (lines roughly 1798..1900). Source the function definitions inline.
ACD_POSTPASS_ENABLED <- TRUE
ACD_POSTPASS_MODE    <- "stratified"
acd_rows <- validation_data$VARIANT == "ACD"
n_acd_local <- sum(acd_rows)
stopifnot(n_acd_local == n_acd)

facs <- fread(factors_path)
facs_wide <- dcast(facs, stratum + level ~ attribute, value.var = "multiplier")

ba_class <- fcase(
  validation_data$BA_t1 < 50,                                "BA1_low",
  validation_data$BA_t1 >= 50  & validation_data$BA_t1 < 100, "BA1_med",
  validation_data$BA_t1 >= 100 & validation_data$BA_t1 < 150, "BA1_high",
  validation_data$BA_t1 >= 150,                              "BA1_vhigh",
  default = NA_character_)

acd_si <- validation_data$SICOND[acd_rows]
si_qs  <- quantile(acd_si, c(1/3, 2/3), na.rm = TRUE)
si_class <- fcase(
  validation_data$SICOND <  si_qs[1],                           "SI_low",
  validation_data$SICOND >= si_qs[1] & validation_data$SICOND < si_qs[2], "SI_med",
  validation_data$SICOND >= si_qs[2],                           "SI_high",
  default = NA_character_)

ft_group <- fcase(
  validation_data$FORTYPCD %in% c(101:109),  "white-red-jack-pine",
  validation_data$FORTYPCD %in% c(121:129),  "spruce-fir",
  validation_data$FORTYPCD %in% c(401:409),  "oak-hickory",
  validation_data$FORTYPCD %in% c(501:509),  "oak-pine",
  validation_data$FORTYPCD %in% c(801:809),  "maple-beech-birch",
  validation_data$FORTYPCD %in% c(901:909),  "aspen-birch",
  default = "nonstocked-other")

int_class <- as.character(pmin(pmax(round(validation_data$interval_years), 4L), 7L))

ft_tab    <- facs_wide[stratum == "FT_GROUP"]
ba_tab    <- facs_wide[stratum == "BA_t1_class"]
si_tab    <- facs_wide[stratum == "SI_tercile"]
int_tab   <- facs_wide[stratum == "interval_years"]

get_mult <- function(tab, levels, attr) {
  idx <- match(levels, tab$level)
  out <- tab[[attr]][idx]
  out[is.na(out)] <- 1.0
  out
}

gm4 <- function(a, b, c, d) (a * b * c * d) ^ 0.25
mult_BAPH  <- gm4(get_mult(ft_tab,  ft_group, "BAPH"),
                  get_mult(ba_tab,  ba_class, "BAPH"),
                  get_mult(si_tab,  si_class, "BAPH"),
                  get_mult(int_tab, int_class,"BAPH"))
mult_TPA   <- gm4(get_mult(ft_tab,  ft_group, "TPA"),
                  get_mult(ba_tab,  ba_class, "TPA"),
                  get_mult(si_tab,  si_class, "TPA"),
                  get_mult(int_tab, int_class,"TPA"))
mult_QMD   <- gm4(get_mult(ft_tab,  ft_group, "QMD"),
                  get_mult(ba_tab,  ba_class, "QMD"),
                  get_mult(si_tab,  si_class, "QMD"),
                  get_mult(int_tab, int_class,"QMD"))
mult_TOPHT <- gm4(get_mult(ft_tab,  ft_group, "TOPHT"),
                  get_mult(ba_tab,  ba_class, "TOPHT"),
                  get_mult(si_tab,  si_class, "TOPHT"),
                  get_mult(int_tab, int_class,"TOPHT"))

cat("Stratum classification (ACD rows):\n")
print(data.table(
  BA_t1      = validation_data$BA_t1[acd_rows],
  SICOND     = validation_data$SICOND[acd_rows],
  FORTYPCD   = validation_data$FORTYPCD[acd_rows],
  interval   = validation_data$interval_years[acd_rows],
  ba_class   = ba_class[acd_rows],
  si_class   = si_class[acd_rows],
  ft_group   = ft_group[acd_rows],
  int_class  = int_class[acd_rows],
  m_BAPH     = round(mult_BAPH[acd_rows],  4),
  m_TPA      = round(mult_TPA[acd_rows],   4),
  m_QMD      = round(mult_QMD[acd_rows],   4),
  m_TOPHT    = round(mult_TOPHT[acd_rows], 4)
))

validation_data[acd_rows, `:=`(
  BA_pred_calib  = BA_pred_calib  * mult_BAPH[acd_rows],
  TPA_pred_calib = TPA_pred_calib * mult_TPA[acd_rows],
  QMD_pred_calib = QMD_pred_calib * mult_QMD[acd_rows]
)]
if ("HT_top_calib" %in% names(validation_data)) {
  validation_data[acd_rows, HT_top_calib := HT_top_calib * mult_TOPHT[acd_rows]]
}

cat("\nAFTER post-pass:\n")
print(validation_data[, .(VARIANT, BA_pred_calib, TPA_pred_calib, QMD_pred_calib, HT_top_calib)])

cat("\nChecks:\n")
ne_unchanged <- all(
  validation_data[VARIANT == "NE", BA_pred_calib]  == 100.0 &
  validation_data[VARIANT == "NE", TPA_pred_calib] == 300.0 &
  validation_data[VARIANT == "NE", QMD_pred_calib] == 7.0 &
  validation_data[VARIANT == "NE", HT_top_calib]   == 60.0)
cat("- NE rows untouched:    ", ne_unchanged, "\n")

acd_changed <- all(validation_data[VARIANT == "ACD", BA_pred_calib] != 100.0)
cat("- ACD BA all modified:  ", acd_changed, "\n")

no_NAs <- !anyNA(validation_data[, .(BA_pred_calib, TPA_pred_calib, QMD_pred_calib, HT_top_calib)])
cat("- No NAs introduced:    ", no_NAs, "\n")

multipliers_finite <- all(is.finite(mult_BAPH[acd_rows]) & is.finite(mult_TPA[acd_rows]) &
                          is.finite(mult_QMD[acd_rows]) & is.finite(mult_TOPHT[acd_rows]))
cat("- All multipliers finite:", multipliers_finite, "\n")

cat("\nSummary multipliers (mean across ACD rows):\n")
cat(sprintf("  BAPH=%.4f  TPA=%.4f  QMD=%.4f  TOPHT=%.4f\n",
  mean(mult_BAPH[acd_rows]),  mean(mult_TPA[acd_rows]),
  mean(mult_QMD[acd_rows]),   mean(mult_TOPHT[acd_rows])))

if (ne_unchanged && acd_changed && no_NAs && multipliers_finite) {
  cat("\nSMOKE TEST: PASS\n")
  quit(status = 0)
} else {
  cat("\nSMOKE TEST: FAIL\n")
  quit(status = 1)
}
