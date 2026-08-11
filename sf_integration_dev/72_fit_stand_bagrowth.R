#!/usr/bin/env Rscript
# =============================================================================
# 72_fit_stand_bagrowth.R  --  STAND-LEVEL basal-area growth model for the CONUS
# arms (symmetric to 71_fit_stand_survival.R). Provides the stand BA-increment
# TARGET that reconciles (disaggregates) the summed tree diameter growth via
# stand_constraint.py::stand_disaggregate_bagrowth (a proportional dg scale).
#
#   response  = stand PERIODIC ANNUAL BA INCREMENT (metric, m2/ha/yr), modeled on
#               the log scale: ln(bai_stand) ~ Normal (a strictly-positive,
#               right-skewed quantity; log keeps predictions positive).
#   eta_bai   = a0 + a_rd*rd + a_lnqmd*ln_qmd + a_ba*ba_metric
#               + a_bgi*bgi + a_bgi2*max(bgi-knot,0)
#               + a_trt*trt_decay + a_dstrb*dstrb_decay + (1 | EPA_L1)
#   BA-growth target (per interval) = exp(eta_bai) * YEARS
#
# Same stand-remeasurement key and covariate construction as 71 (rd = additive
# SDI / Bayesian SDImax_brms, ln_qmd, ba_metric, bgi 2-piece spline, mgmt/dstrb
# decays), hierarchical by EPA_L1, brms, seeded 20260702, fault-tolerant, gc().
# Writes stand_bagrowth_bundle.json.
#
# Usage:
#   Rscript 72_fit_stand_bagrowth.R --n_sub=120000 \
#     --pairs=<rds> --out_dir=<dir> --tau_m=10 --tau_d=15
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
TAU_M <- as.numeric(ga("tau_m", "10"))     # [ASSUMPTION: mgmt e-folding ~10 yr, mirrors 70/71]
TAU_D <- as.numeric(ga("tau_d", "15"))     # [ASSUMPTION: dstrb e-folding ~15 yr, mirrors 70/71]
PAIRS <- ga("pairs", "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR <- ga("out_dir", "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_bagrowth")
N_SUB <- as.integer(ga("n_sub", "120000"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 72_fit_stand_bagrowth.R ==\n tau_m:", TAU_M, " tau_d:", TAU_D, " n_sub:", N_SUB, "\n")

d <- trycatch_run(as.data.table(readRDS(PAIRS)), "read pairs")
if (is.null(d)) quit(status = 1)
cat(" pairs rows:", nrow(d), "\n")

# ---- row-level filter: live-at-t1 trees with a valid remeasurement ----------
# BA increment is built from surviving trees' DBH growth; keep trees alive at t1
# with finite DBH1/DBH2 over a sane interval. Trees that died contribute their
# t1 BA to the stand at t1 but no positive growth (their loss is the survival
# model's job); for a stand BA GROWTH (gross of ingrowth, net of mortality on
# survivors) we sum the per-tree BA increment of SURVIVORS (TREESTATUS2==1).
if (!"TREESTATUS1" %in% names(d)) d[, TREESTATUS1 := 1L]
d <- d[TREESTATUS1 == 1 & !is.na(TREESTATUS2) &
       is.finite(TPA1) & TPA1 > 0 & is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
       is.finite(DBH1) & DBH1 >= 2.54]
cat(" live-at-t1 tree rows:", nrow(d), "\n")

# ---- STAND KEY (identical to 71_fit_stand_survival.R) -----------------------
d[, stand_key := paste(PLT_CN_cond1, CONDID_cond1, INVYR1, INVYR2, sep = "_")]

# ---- per-tree BA increment (metric): m2 per stem, survivors only -------------
# BA (m2) = pi/4 * (DBH_cm/100)^2. Survivors (TREESTATUS2==1) with a valid DBH2
# grow; dead/removed contribute 0 growth. TPA-weight to per-ha increment.
K <- pi / 4
d[, ba1_stem := K * (DBH1 / 100)^2]
d[, dbh2_use := ifelse(TREESTATUS2 == 1 & is.finite(DBH2) & DBH2 >= DBH1, DBH2, DBH1)]
d[, ba2_stem := K * (dbh2_use / 100)^2]
d[, ba_incr_stem := pmax(ba2_stem - ba1_stem, 0)]     # non-negative growth

# ---- aggregate tree pairs -> stand-remeasurement records ---------------------
d[, died := as.numeric(TREESTATUS2 == 2)]
agg <- d[, .(
  YEARS      = YEARS[1],
  ba_incr    = sum(TPA1 * ba_incr_stem),              # stand BA increment (m2/ha over interval)
  ntrees     = sum(TPA1),
  ba_metric  = weighted.mean(BA1, TPA1, na.rm = TRUE),
  qmd        = weighted.mean(QMD1, TPA1, na.rm = TRUE),
  sdi_add    = weighted.mean(sdi_additive1, TPA1, na.rm = TRUE),
  sdimax     = weighted.mean(SDImax_brms, TPA1, na.rm = TRUE),
  rd_row     = weighted.mean(rd_add, TPA1, na.rm = TRUE),
  bgi        = weighted.mean(bgi, TPA1, na.rm = TRUE),
  yst        = suppressWarnings(min(years_since_trt,   na.rm = TRUE)),
  ysd        = suppressWarnings(min(years_since_dstrb, na.rm = TRUE)),
  L1         = EPA_L1_CODE[1]
), by = stand_key]
rm(d); gc()
cat(" stand-remeasurement records:", nrow(agg), "\n")

# ---- stand covariates (identical construction to 71) ------------------------
agg[, rd := ifelse(is.finite(sdi_add) & is.finite(sdimax) & sdimax > 0, sdi_add / sdimax, rd_row)]
# periodic annual BA increment (m2/ha/yr), response on the log scale
agg[, bai_annual := ba_incr / YEARS]
agg <- agg[is.finite(YEARS) & is.finite(bai_annual) & bai_annual > 0 &
           is.finite(rd) & rd > 0 & rd < 2 &
           is.finite(qmd) & qmd > 0 & is.finite(ba_metric) & ba_metric >= 0]
agg[, ln_qmd := log(qmd)]
agg[, ln_bai := log(bai_annual)]
bgi_med <- median(agg$bgi[is.finite(agg$bgi)], na.rm = TRUE)
agg[!is.finite(bgi), bgi := bgi_med]
agg[, bgi_b2 := pmax(bgi - bgi_med, 0)]
agg[, trt_active   := as.integer(is.finite(yst) & yst >= 0)]
agg[, dstrb_active := as.integer(is.finite(ysd) & ysd >= 0)]
agg[, trt_decay   := ifelse(trt_active   == 1, exp(-pmin(yst, 100) / TAU_M), 0)]
agg[, dstrb_decay := ifelse(dstrb_active == 1, exp(-pmin(ysd, 100) / TAU_D), 0)]
agg[, L1 := as.character(L1)]
agg <- agg[!is.na(L1) & L1 != ""]
cat(" modeling-ready stand records:", nrow(agg), "\n")
cat(" mean stand BAI (m2/ha/yr):", round(mean(agg$bai_annual), 4),
    " median:", round(median(agg$bai_annual), 4), "\n")

# ---- subsample (oversample event stands, as in 71) --------------------------
ev <- agg[trt_active == 1 | dstrb_active == 1]
bg <- agg[trt_active == 0 & dstrb_active == 0]
n_ev <- min(nrow(ev), as.integer(N_SUB * 0.5))
n_bg <- min(nrow(bg), N_SUB - n_ev)
ds <- rbind(ev[sample(.N, n_ev)], bg[sample(.N, n_bg)])
cat(sprintf(" fit sample: %d (events %d, baseline %d) of %d eligible\n",
            nrow(ds), n_ev, n_bg, nrow(agg)))
rm(agg, ev, bg); gc()

# ---- Bayesian fit (brms): Gaussian on ln(BAI) -------------------------------
ok <- trycatch_run({ suppressPackageStartupMessages(library(brms)); TRUE }, "load brms")
if (is.null(ok)) quit(status = 1)

FIXED <- c("rd", "ln_qmd", "ba_metric", "bgi", "bgi_b2", "trt_decay", "dstrb_decay")
rhs <- paste(FIXED, collapse = " + ")
form <- bf(as.formula(paste0("ln_bai ~ ", rhs, " + (1 | L1)")))
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"),
            set_prior("student_t(3,0,2.5)", class = "sigma"))
fit <- trycatch_run(
  brm(form, data = ds, family = gaussian(), prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, control = list(adapt_delta = 0.9)),
  "brms stand-bagrowth fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
sm <- as.data.frame(fx); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, "stand_bagrowth_fixed.csv"))
sigma_resid <- trycatch_run(unname(summary(fit)$spec_pars["sigma", "Estimate"]), "sigma")

fe_names <- rownames(fx)
fixed_list <- setNames(
  lapply(fe_names, function(nm) list(mean = unname(fx[nm, "Estimate"]),
                                     sd   = unname(fx[nm, "Est.Error"]),
                                     q2.5 = unname(fx[nm, "Q2.5"]),
                                     q97.5= unname(fx[nm, "Q97.5"]))),
  fe_names)
re_table <- if (!is.null(re)) list(level = rownames(re),
                                   mean  = unname(re[, "Estimate"]),
                                   sd    = unname(re[, "Est.Error"])) else NULL
sd_L1 <- trycatch_run({ vc <- summary(fit)$random$L1; unname(vc["sd(Intercept)", "Estimate"]) }, "sd L1")

bundle <- list(
  model      = "stand-level basal-area growth (log periodic annual BA increment)",
  scale      = "Gaussian on ln(BAI); BAI = stand basal-area increment m2/ha/yr; predict exp(eta) then * YEARS for the interval target",
  response   = "ln(bai_annual); bai_annual = sum_i TPA_i * max(BA2_i - BA1_i, 0) / YEARS (survivors)",
  stand_key  = "PLT_CN_cond1 x CONDID_cond1 x INVYR1 x INVYR2 (matches 71_fit_stand_survival.R)",
  target     = "stand BA increment over interval = exp(linpred) * YEARS; feeds stand_constraint.py::stand_disaggregate_bagrowth (proportional tree-dg scale)",
  covariates = list(
    rd          = "stand SDI / SDIMAX (sdi_additive1 / SDImax_brms)",
    ln_qmd      = "log stand QMD (cm)",
    ba_metric   = "stand basal area (m2/ha)",
    bgi         = "stand mean bgi; 2-piece spline knot at median",
    bgi_knot    = bgi_med,
    trt_decay   = "exp(-min(years_since_trt)/tau_m)",
    dstrb_decay = "exp(-min(years_since_dstrb)/tau_d)"),
  tau_m = TAU_M, tau_d = TAU_D,
  n_stand_records = nrow(ds),
  fixed_effects = fixed_list,
  sigma_resid = sigma_resid,
  sd_L1 = sd_L1,
  re_L1 = re_table,
  seed = 20260702,
  notes = "Symmetric to 71_fit_stand_survival.R. Stand BA-growth target reconciles the summed tree DG via a proportional factor gamma on tree dg (BA increment is additive)."
)
write_json(bundle, file.path(OUT_DIR, "stand_bagrowth_bundle.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")

cat("\nStand-BA-growth fixed effects (ln BAI scale):\n")
print(round(fx[, c("Estimate", "Est.Error", "Q2.5", "Q97.5"), drop = FALSE], 4))
cat("\nsigma:", round(ifelse(is.null(sigma_resid), NA, sigma_resid), 4),
    " sd(1|L1):", round(ifelse(is.null(sd_L1), NA, sd_L1), 4), "\n")
gc(); cat("\nDONE_STAND_BAGROWTH\n")
