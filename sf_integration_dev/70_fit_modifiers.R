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
# Each component defines: y (response on its transformed/link scale), the row
# filters, the base-control columns, BASE_CTRL (their names), and FAMILY.
# The MODIFIER design (trt_decay, dstrb_decay, bgi, bgi_b2, (1|L1)) and the
# bundle output structure are IDENTICAL across components; only the response +
# base controls below differ. Response definitions + covariate forms mirror
# ~/fvs-conus/dev_sf_integration/benchmark_sf_vs_legA.R (prep_hg / prep_htdbh /
# prep_hcb / prep_cr2) and benchmark_mort_legA.R (survival Bernoulli form).
logit <- function(p) log(p / (1 - p))
# FAMILY is chosen here as a NAME string ("gaussian"/"bernoulli") and the brms
# family object is constructed after library(brms) loads (bernoulli() is a brms
# function, not base R, so it cannot be called before the package is attached).
FAMILY <- "gaussian"          # overridden to "bernoulli" for mort below
if (COMPONENT == "dg") {
  # diameter growth: log(annual dDBH); base controls ln_dbh/ln_cr_adj/comp
  d[, resp := (DBH2 - DBH1) / YEARS]
  d[, y := log(pmax(resp, 1e-4))]
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & is.finite(resp) & resp > 0.01 & resp < 5.0]
  d[, ln_dbh := log(DBH1)]
  d[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  d[, comp := (BAL_SW1 + BAL_HW1)]
  BASE_CTRL <- c("ln_dbh", "ln_cr_adj", "comp")

} else if (COMPONENT == "hg") {
  # height growth: log(annual dHT) (gaussian); prep_hg covariate forms
  d[, resp := (HT2 - HT1) / YEARS]
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & is.finite(HT1) & HT1 > 1.5 &
         is.finite(HT2) & HT2 > 1.5 & is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
         is.finite(resp) & resp > 0.001 & resp < 5.0]
  d[, y := log(resp)]
  d[, ln_dbh := log(DBH1)]
  d[, ln_ht := log(pmax(HT1, 1.5))]
  d[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  d[, bal_log := log((BAL_SW1 + BAL_HW1) + 5)]
  d[, ba_metric := BA1 * 0.2296]
  BASE_CTRL <- c("ln_dbh", "ln_ht", "ln_cr_adj", "bal_log", "ba_metric")

} else if (COMPONENT == "htdbh") {
  # height-diameter: log(HT-1.37) (gaussian); prep_htdbh covariate forms
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & DBH1 < 250 & is.finite(HT1) & HT1 >= 1.37 & HT1 < 85 &
         is.finite(BA1) & BA1 > 0 & is.finite(BAL_SW1) & is.finite(BAL_HW1) &
         is.finite(sdi_additive1) & is.finite(SDImax_brms) & SDImax_brms > 0]
  d[, resp := HT1]
  d[, y := log(pmax(HT1 - 1.37, 0.01))]
  d[, bal := (BAL_SW1 + BAL_HW1)]
  d[, sqrt_ba := sqrt(BA1 * 0.2296)]
  d[, rd_ratio := sdi_additive1 / SDImax_brms]
  d[, ba_metric := BA1 * 0.2296]
  d[, ba_x_rd := ba_metric * rd_ratio]
  d[, bal_x_rd := bal * rd_ratio]
  d[, inv_dbh := 1 / (DBH1 + 1)]
  d <- d[is.finite(rd_ratio) & rd_ratio > 0 & rd_ratio < 2]
  BASE_CTRL <- c("bal", "sqrt_ba", "ba_x_rd", "bal_x_rd", "inv_dbh")

} else if (COMPONENT == "hcb") {
  # height-crown-base: logit((1-CR1) mapped to (0,1)) (gaussian on logit scale);
  # prep_hcb response = 1-CR1, covariate forms from prep_hcb
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & is.finite(HT1) & HT1 > 1 &
         is.finite(CR1) & CR1 > 0 & CR1 < 1 & is.finite(BA1) & BA1 >= 0 &
         is.finite(BAL_SW1) & is.finite(BAL_HW1)]
  d[, resp := 1 - CR1]
  d <- d[resp > 0.01 & resp < 0.99]
  d[, y := logit(pmin(pmax((resp - 0.001) / 0.998, 1e-6), 1 - 1e-6))]
  d[, ln_ht := log(pmax(HT1, 1.5))]
  d[, ln_dbh := log(DBH1)]
  d[, bal_over_ht := (BAL_SW1 + BAL_HW1) / (HT1 + 1)]
  d[, sqrt_ba := sqrt(BA1 * 0.2296)]
  BASE_CTRL <- c("ln_ht", "ln_dbh", "bal_over_ht", "sqrt_ba")

} else if (COMPONENT == "cr") {
  # crown recession: CR2-direct logit(CR2) (gaussian on logit scale); prep_cr2
  # covariate forms (logit(CR1) + DBH + DBH^2 + BA + BAL controls)
  d <- d[is.finite(DBH1) & DBH1 >= 2.54 & is.finite(CR1) & CR1 > 0 & CR1 < 1 &
         is.finite(CR2) & CR2 > 0 & CR2 < 1 & is.finite(BA1) & BA1 >= 0 &
         is.finite(BAL_SW1) & is.finite(BAL_HW1)]
  d[, resp := CR2]
  d[, y := logit(pmin(pmax(CR2, 1e-4), 1 - 1e-4))]
  d[, cr1_logit := logit(pmin(pmax(CR1, 1e-4), 1 - 1e-4))]
  d[, dbh := DBH1]
  d[, dbh_sq := DBH1^2]
  d[, ba_metric := BA1 * 0.2296]
  d[, bal_metric := (BAL_SW1 + BAL_HW1)]
  BASE_CTRL <- c("cr1_logit", "dbh", "dbh_sq", "ba_metric", "bal_metric")

} else if (COMPONENT == "mort") {
  # mortality: continuous-time exponential-hazard survival, CONSISTENT with the
  # base model survival_unified_v2_crz (log S = -exp(-eta) * T_years). That maps
  # to a cloglog GLM of the mortality EVENT with a log(YEARS) exposure offset:
  #   P(death) = 1 - exp(-exp(eta_h) * T),  cloglog(P_death) = eta_h + log(T).
  # The modifier terms enter the log-hazard eta_h and multiply the base hazard by
  # exp(mod_eta), matching how the base survival model applies effects.
  FAMILY <- "bernoulli_cloglog"
  if (!"TREESTATUS1" %in% names(d)) d[, TREESTATUS1 := 1L]
  d <- d[TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1, 2) &
         is.finite(DBH1) & DBH1 >= 2.54 & is.finite(CR1) & CR1 > 0 & CR1 <= 1 &
         is.finite(BA1) & BA1 >= 0 & is.finite(BAL_SW1) & is.finite(BAL_HW1) &
         is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
         is.finite(sdi_additive1) & is.finite(SDImax_brms) & SDImax_brms > 0]
  d[, y := as.integer(TREESTATUS2 == 2)]      # response: mortality event (died)
  d[, resp := y]
  d[, dbh := DBH1]
  d[, dbh_sq := DBH1^2]
  cr_m <- mean(d$CR1, na.rm = TRUE); cr_s <- sd(d$CR1, na.rm = TRUE); if (!is.finite(cr_s) || cr_s == 0) cr_s <- 1
  d[, cr_z := (CR1 - cr_m) / cr_s]
  d[, cr_z_sq := cr_z^2]
  d[, rd_ratio := sdi_additive1 / SDImax_brms]
  d[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]
  d[, bal_metric := (BAL_SW1 + BAL_HW1)]
  # exposure offset with a nominal baseline log-hazard folded in (~1.8%/yr) so
  # the intercept starts near a realistic mortality rate; without this the
  # cloglog inverse of log(YEARS) is ~1 at init and survivors give log(0)=-Inf.
  d[, log_years := log(YEARS) - 3.9]          # log(YEARS) + log(~0.02) baseline hazard
  d <- d[is.finite(sqrt_ba_rd) & is.finite(log_years)]
  BASE_CTRL <- c("dbh", "dbh_sq", "cr_z", "cr_z_sq", "bal_metric", "sqrt_ba_rd")

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
d <- d[!is.na(L1) & L1 != ""]
# drop rows with any non-finite base control (component-specific set)
for (cc in BASE_CTRL) d <- d[is.finite(get(cc))]

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

# base controls (component-specific) + IDENTICAL modifier terms + L1 RE.
rhs <- paste(c(BASE_CTRL,                            # base controls (partial out size/competition)
               "trt_decay", "dstrb_decay",          # management + disturbance modifiers
               "bgi", "bgi_b2"),                    # climate driver modifier
             collapse = " + ")
off <- if (identical(FAMILY, "bernoulli_cloglog")) " + offset(log_years)" else ""
form <- bf(as.formula(paste0("y ~ ", rhs, " + (1 | L1)", off)))   # + ecoregion RE (+ exposure offset for mort)
is_bern <- FAMILY %in% c("bernoulli", "bernoulli_cloglog")
# construct the brms family object now that the package is attached; mortality
# uses a cloglog link with a log(YEARS) exposure offset -> modifier on log-hazard
fam_obj <- if (identical(FAMILY, "bernoulli_cloglog")) bernoulli(link = "cloglog") else if (is_bern) bernoulli() else gaussian()
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"))
if (!is_bern) priors <- c(priors, set_prior("student_t(3,0,1)", class = "sigma"))
fit <- trycatch_run(
  brm(form, data = ds, family = fam_obj, prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, init = 0, control = list(adapt_delta = 0.9)),
  "brms fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)                    # posterior mean + CI for fixed effects
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
sm <- as.data.frame(fx); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, paste0(COMPONENT, "_modifier_fixed.csv")))

# modifier bundle: only the modifier terms (exclude base controls + global Intercept)
mod_terms <- c("trt_decay", "dstrb_decay", "bgi", "bgi_b2")
# residual scale: gaussian components carry a sigma spec_par; the bernoulli
# (mortality) fit has none, so tag sigma_resid = NA on the response (logit) scale.
sigma_resid <- if (is_bern) NA_real_ else {
  sp <- summary(fit)$spec_pars
  if (!is.null(sp) && "sigma" %in% rownames(sp)) as.numeric(sp["sigma", "Estimate"]) else NA_real_
}
bundle <- list(
  component = COMPONENT, form = "multiplicative log-modifier",
  family = FAMILY,
  response_scale = if (identical(FAMILY, "bernoulli_cloglog")) "log-hazard (cloglog, log-years exposure offset); multiplier = exp(mod_eta) on the base hazard, consistent with survival_unified_v2_crz" else if (is_bern) "logit(survival)" else "component transform",
  tau_m = TAU_M, tau_d = TAU_D, bgi_knot = bgi_med,
  base_controls = BASE_CTRL,
  sigma_resid = sigma_resid,
  management = list(param = "trt_decay", coef = unname(fx["trt_decay", "Estimate"]),
                    coef_sd = unname(fx["trt_decay", "Est.Error"]),
                    decay = "exp(-years_since_trt / tau_m)"),
  disturbance = list(param = "dstrb_decay", coef = unname(fx["dstrb_decay", "Estimate"]),
                     coef_sd = unname(fx["dstrb_decay", "Est.Error"]),
                     decay = "exp(-years_since_dstrb / tau_d)"),
  driver_bgi = list(b1 = unname(fx["bgi", "Estimate"]), b1_sd = unname(fx["bgi", "Est.Error"]),
                    b2 = unname(fx["bgi_b2", "Estimate"]), b2_sd = unname(fx["bgi_b2", "Est.Error"]),
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
