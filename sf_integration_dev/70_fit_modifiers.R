#!/usr/bin/env Rscript
# =============================================================================
# 70_fit_modifiers.R  --  Bayesian modifier layer for the three CONUS arms.
#
# Fits a shared MULTIPLICATIVE modifier on log-growth that any base arm
# (conus_organon / conus / conus_sf) applies on top of its native prediction:
#
#   growth = base_arm_prediction * exp(mod_eta)
#   mod_eta = a_trt * trt_decay(years_since_trt)          # management
#           + a_dstrb * dstrb_decay(years_since_dstrb)     # disturbance
#           + b_bgi1 * bgi + b_bgi2 * (bgi - knot)_+       # climate driver
#           + z_L1[ecoregion]                              # hierarchical RE
#
# Baseline (no recent event) -> trt_decay = dstrb_decay = 0 -> multiplier 1.
# Estimated as the PARTIAL effect of management / disturbance / drivers holding
# base size and competition controls fixed, so the same modifier layer is valid
# regardless of which base arm is selected.
#
# Usage:
#   Rscript 70_fit_modifiers.R --component=dg --tau_m=10 --tau_d=15 \
#     --pairs=<rds> --out_dir=<dir> --n_sub=200000
#
# Author: A. Weiskittel + Claude (OODA autopilot)  Date: 2026-07-02
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
set.seed(20260702)
ELOG <- "error_log.txt"
elog <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time()), msg), file = ELOG, append = TRUE)
trycatch_run <- function(expr, what) tryCatch(expr, error = function(e) { elog(paste(what, ":", conditionMessage(e))); NULL })

args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d = NULL) { m <- grep(paste0("^--", n, "="), args, value = TRUE); if (length(m)) sub(paste0("^--", n, "="), "", m[1]) else d }
COMPONENT <- ga("component", "dg")
TAU_M <- as.numeric(ga("tau_m", "10"))     # [ASSUMPTION: management effect e-folding ~10 yr]
TAU_D <- as.numeric(ga("tau_d", "15"))     # [ASSUMPTION: disturbance effect e-folding ~15 yr]
PAIRS <- ga("pairs", "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR <- ga("out_dir", "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/modifiers")
N_SUB <- as.integer(ga("n_sub", "200000"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 70_fit_modifiers.R ==\n component:", COMPONENT, " tau_m:", TAU_M, " tau_d:", TAU_D, "\n")

d <- trycatch_run(as.data.table(readRDS(PAIRS)), "read pairs")
if (is.null(d)) quit(status = 1)

# ---- response + base controls (component-specific) --------------------------
if (COMPONENT == "dg") {
  d[, resp := (DBH2 - DBH1) / YEARS]
  d[, y := log(pmax(resp, 1e-4))]
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & is.finite(resp) & resp > 0.01 & resp < 5.0]
  d[, ln_dbh := log(DBH1)]
  d[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  d[, comp := (BAL_SW1 + BAL_HW1)]
} else {
  elog(paste("component not yet parameterized:", COMPONENT)); quit(status = 1)
}

# ---- modifier design: decay basis for management + disturbance --------------
d[, trt_active   := as.integer(is.finite(years_since_trt)   & years_since_trt   >= 0)]
d[, dstrb_active := as.integer(is.finite(years_since_dstrb) & years_since_dstrb >= 0)]
d[, trt_decay   := ifelse(trt_active   == 1, exp(-pmin(years_since_trt,   100) / TAU_M), 0)]
d[, dstrb_decay := ifelse(dstrb_active == 1, exp(-pmin(years_since_dstrb, 100) / TAU_D), 0)]
# climate driver spline (2-piece), knot at median bgi among finite
bgi_med <- median(d$bgi[is.finite(d$bgi)], na.rm = TRUE)
d[!is.finite(bgi), bgi := bgi_med]
d[, bgi_b2 := pmax(bgi - bgi_med, 0)]
d[, L1 := as.character(EPA_L1_CODE)]
d <- d[!is.na(L1) & L1 != "" & is.finite(ln_dbh) & is.finite(ln_cr_adj) & is.finite(comp)]

# stratified subsample: oversample event pairs so the modifier signal is seen
ev <- d[trt_active == 1 | dstrb_active == 1]
bg <- d[trt_active == 0 & dstrb_active == 0]
n_ev <- min(nrow(ev), as.integer(N_SUB * 0.5))
n_bg <- min(nrow(bg), N_SUB - n_ev)
ds <- rbind(ev[sample(.N, n_ev)], bg[sample(.N, n_bg)])
cat(sprintf(" fit sample: %d (events %d, baseline %d) of %d eligible\n", nrow(ds), n_ev, n_bg, nrow(d)))
rm(d, ev, bg); gc()

# ---- Bayesian fit (brms) : base controls + modifier terms + L1 RE ----------
ok <- trycatch_run({ suppressPackageStartupMessages(library(brms)); TRUE }, "load brms")
if (is.null(ok)) quit(status = 1)

form <- bf(y ~ ln_dbh + ln_cr_adj + comp +          # base controls (partial out size/competition)
             trt_decay + dstrb_decay +              # management + disturbance modifiers
             bgi + bgi_b2 +                         # climate driver modifier
             (1 | L1))                              # ecoregion RE
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"),
            set_prior("student_t(3,0,1)", class = "sigma"))
fit <- trycatch_run(
  brm(form, data = ds, family = gaussian(), prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, control = list(adapt_delta = 0.9)),
  "brms fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)                    # posterior mean + CI for fixed effects
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
sm <- as.data.frame(fx); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, paste0(COMPONENT, "_modifier_fixed.csv")))

# modifier bundle: only the modifier terms (exclude base controls + global Intercept)
mod_terms <- c("trt_decay", "dstrb_decay", "bgi", "bgi_b2")
bundle <- list(
  component = COMPONENT, form = "multiplicative log-modifier",
  tau_m = TAU_M, tau_d = TAU_D, bgi_knot = bgi_med,
  management = list(param = "trt_decay", coef = unname(fx["trt_decay", "Estimate"]),
                    decay = "exp(-years_since_trt / tau_m)"),
  disturbance = list(param = "dstrb_decay", coef = unname(fx["dstrb_decay", "Estimate"]),
                     decay = "exp(-years_since_dstrb / tau_d)"),
  driver_bgi = list(b1 = unname(fx["bgi", "Estimate"]), b2 = unname(fx["bgi_b2", "Estimate"]),
                    knot = bgi_med),
  re_L1 = if (!is.null(re)) list(level = rownames(re), mean = unname(re[, "Estimate"])) else NULL,
  notes = "Multiplier = exp(mod_eta). mod_eta = coef_mgmt*trt_decay + coef_dstrb*dstrb_decay + b1*bgi + b2*(bgi-knot)_+ + z_L1. Baseline (no event) => trt/dstrb decay 0."
)
write_json(bundle, file.path(OUT_DIR, paste0(COMPONENT, "_modifier_bundle.json")),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
cat("\nModifier coefficients (log scale):\n")
print(round(fx[mod_terms, c("Estimate", "Q2.5", "Q97.5"), drop = FALSE], 4))
cat("\nMANAGEMENT multiplier at years_since_trt=0:", round(exp(fx["trt_decay","Estimate"]), 3), "\n")
cat("DISTURBANCE multiplier at years_since_dstrb=0:", round(exp(fx["dstrb_decay","Estimate"]), 3), "\n")
gc(); cat("\nDONE_MODFIT\n")
