#!/usr/bin/env Rscript
# =============================================================================
# 71_fit_stand_survival.R  --  STAND-LEVEL survival model for the CONUS arms.
#
# Fits a stand-level continuous-time exponential-hazard survival model that is
# CONSISTENT with the tree model's scale (survival_unified_v2_crz_dcch.stan and
# the mortality modifier in 70_fit_modifiers.R):
#
#   S_stand(T)      = exp(-H_stand * T_years)
#   H_stand         = exp(-eta_stand)                 # annual stand hazard
#   cloglog(M_stand)= log(H_stand) + log(T) = -eta_stand + log(T)   (exposure offset)
#   eta_stand       = f(rd, ln_qmd, ba_metric, bgi, mgmt/dstrb decays) + (1 | EPA_L1)
#
# Bigger eta_stand -> lower hazard -> higher stand survival. This is the SAME
# cloglog + log(YEARS) exposure form the tree model uses, so the stand target
# M_stand = 1 - S_stand(T) can be used to reconcile (disaggregate) the summed
# tree hazards via the proportional-hazard kappa solve in stand_constraint.py.
#
# The stand records are built by aggregating the tree remeasurement pairs to a
# stand-remeasurement key (see STAND KEY below). Each stand-remeasurement is a
# binomial trial: trials = TPA-weighted live trees at time1, deaths = TPA-weighted
# count of trees that died (TREESTATUS2 == 2) over the interval.
#
# Usage:
#   Rscript 71_fit_stand_survival.R --n_sub=120000 \
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
TAU_M <- as.numeric(ga("tau_m", "10"))     # [ASSUMPTION: management effect e-folding ~10 yr, mirrors 70_fit_modifiers.R]
TAU_D <- as.numeric(ga("tau_d", "15"))     # [ASSUMPTION: disturbance effect e-folding ~15 yr, mirrors 70_fit_modifiers.R]
PAIRS <- ga("pairs", "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR <- ga("out_dir", "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_survival")
N_SUB <- as.integer(ga("n_sub", "120000"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 71_fit_stand_survival.R ==\n tau_m:", TAU_M, " tau_d:", TAU_D, " n_sub:", N_SUB, "\n")

d <- trycatch_run(as.data.table(readRDS(PAIRS)), "read pairs")
if (is.null(d)) quit(status = 1)
cat(" pairs rows:", nrow(d), "\n")

# ---- row-level filter (same live-at-t1 basis as the tree mortality fit) ------
# TREESTATUS2 == 1 survived, == 2 died (mortality), == 3 removed/cut (excluded,
# not a mortality event). Keep live trees at t1 with a valid t2 live/dead status.
if (!"TREESTATUS1" %in% names(d)) d[, TREESTATUS1 := 1L]
d <- d[TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1, 2) &
       is.finite(TPA1) & TPA1 > 0 & is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
       is.finite(DBH1) & DBH1 >= 2.54]
cat(" live-at-t1 tree rows:", nrow(d), "\n")

# ---- STAND KEY --------------------------------------------------------------
# FIA plot-condition-remeasurement: PLT_CN_cond1 x CONDID_cond1 x INVYR1 x INVYR2.
# Verified: 325,486 unique keys, each with a SINGLE YEARS value (a well-formed
# stand-remeasurement). plot_key alone collapses multiple remeasurements of the
# same plot, so the interval (INVYR1,INVYR2) is included. No fallback needed.
d[, stand_key := paste(PLT_CN_cond1, CONDID_cond1, INVYR1, INVYR2, sep = "_")]

# ---- aggregate tree pairs -> stand-remeasurement records ---------------------
# stand trials  = sum of per-tree TPA at t1 (TPA-weighted live count)
# stand deaths  = TPA-weighted count of trees with TREESTATUS2 == 2
# stand survival is over the interval YEARS (constant within a stand key).
# Stand covariates are TPA-weighted means of the per-row (already stand-level or
# tree-level) fields; rd / SDImax / bgi / BA / QMD are stand-constant, so the
# weighted mean simply returns the stand value. mgmt/dstrb use the MOST RECENT
# event in the stand (min years_since_*), i.e. the strongest active effect.
d[, died := as.numeric(TREESTATUS2 == 2)]
agg <- d[, .(
  YEARS      = YEARS[1],
  trials     = sum(TPA1),
  deaths     = sum(TPA1 * died),
  # stand structure (metric): BA m2/ha, QMD cm, SDImax, relative density rd
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

# ---- stand covariates -------------------------------------------------------
# relative density rd = stand SDI / SDIMAX (Bayesian SDImax_brms). Prefer the
# recomputed rd = sdi_add/sdimax; fall back to the row-level rd_add mean.
agg[, rd := ifelse(is.finite(sdi_add) & is.finite(sdimax) & sdimax > 0, sdi_add / sdimax, rd_row)]
agg <- agg[is.finite(YEARS) & is.finite(trials) & trials > 0 & is.finite(deaths) &
           deaths >= 0 & deaths <= trials &
           is.finite(rd) & rd > 0 & rd < 2 &
           is.finite(qmd) & qmd > 0 & is.finite(ba_metric) & ba_metric >= 0]
agg[, ln_qmd := log(qmd)]
# bgi: fill non-finite with median, 2-piece spline (knot at median) as in 70_*
bgi_med <- median(agg$bgi[is.finite(agg$bgi)], na.rm = TRUE)
agg[!is.finite(bgi), bgi := bgi_med]
agg[, bgi_b2 := pmax(bgi - bgi_med, 0)]
# management / disturbance decays (most-recent event in the stand)
agg[, trt_active   := as.integer(is.finite(yst) & yst >= 0)]
agg[, dstrb_active := as.integer(is.finite(ysd) & ysd >= 0)]
agg[, trt_decay   := ifelse(trt_active   == 1, exp(-pmin(yst, 100) / TAU_M), 0)]
agg[, dstrb_decay := ifelse(dstrb_active == 1, exp(-pmin(ysd, 100) / TAU_D), 0)]
agg[, log_years := log(YEARS) - 3.9]  # + nominal baseline log-hazard (~1.8%/yr) so cloglog init is not P=1                       # exposure offset (cloglog)
agg[, L1 := as.character(L1)]
agg <- agg[!is.na(L1) & L1 != "" & is.finite(log_years)]
# integer trial/death counts for the binomial response
agg[, trials_i := as.integer(round(pmax(trials, 1)))]
agg[, deaths_i := as.integer(pmin(round(deaths), trials_i))]
cat(" modeling-ready stand records:", nrow(agg), "\n")
cat(" mean stand mortality fraction:", round(mean(agg$deaths_i / agg$trials_i), 4), "\n")

# ---- subsample for the Bayesian fit -----------------------------------------
# oversample stands with a recent event so the mgmt/dstrb signal is represented
ev <- agg[trt_active == 1 | dstrb_active == 1]
bg <- agg[trt_active == 0 & dstrb_active == 0]
n_ev <- min(nrow(ev), as.integer(N_SUB * 0.5))
n_bg <- min(nrow(bg), N_SUB - n_ev)
ds <- rbind(ev[sample(.N, n_ev)], bg[sample(.N, n_bg)])
cat(sprintf(" fit sample: %d (events %d, baseline %d) of %d eligible\n",
            nrow(ds), n_ev, n_bg, nrow(agg)))
rm(agg, ev, bg); gc()

# ---- Bayesian fit (brms): cloglog binomial + log(YEARS) exposure offset ------
ok <- trycatch_run({ suppressPackageStartupMessages(library(brms)); TRUE }, "load brms")
if (is.null(ok)) quit(status = 1)

# eta_stand fixed effects: relative density, ln(QMD), BA, bgi spline, mgmt/dstrb
# decays. cloglog(deaths/trials) = eta_h + log(YEARS); eta_h = -eta_stand, i.e.
# the linear predictor here IS the log-hazard. Reported so the disaggregation
# reconciler can reconstruct H_stand = exp(linpred - <structure only>) as needed.
FIXED <- c("rd", "ln_qmd", "ba_metric", "bgi", "bgi_b2", "trt_decay", "dstrb_decay")
rhs <- paste(FIXED, collapse = " + ")
form <- bf(as.formula(paste0("deaths_i | trials(trials_i) ~ ", rhs,
                             " + (1 | L1) + offset(log_years)")))
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"),
            set_prior("normal(0,2)", class = "Intercept"))
fit <- trycatch_run(
  brm(form, data = ds, family = binomial(link = "cloglog"), prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, control = list(adapt_delta = 0.9)),
  "brms stand-survival fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)                                     # posterior mean + SD + CI
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
sm <- as.data.frame(fx); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, "stand_survival_fixed.csv"))

# ---- bundle: fixed-effect posterior means + SDs + RE table -------------------
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
  model      = "stand-level survival (continuous-time exponential hazard)",
  scale      = "cloglog binomial with log(YEARS) exposure offset; linear predictor = log-hazard, consistent with survival_unified_v2_crz and the mortality modifier (70_fit_modifiers.R)",
  response   = "deaths_i | trials(trials_i); trials = TPA-weighted live trees at t1, deaths = TPA-weighted TREESTATUS2==2",
  stand_key  = "PLT_CN_cond1 x CONDID_cond1 x INVYR1 x INVYR2 (FIA plot-condition-remeasurement)",
  hazard     = "H_stand = exp(linpred); S_stand(T) = exp(-H_stand * T); M_stand = 1 - S_stand(T)",
  covariates = list(
    rd          = "stand SDI / SDIMAX (sdi_additive1 / SDImax_brms, Bayesian SDImax)",
    ln_qmd      = "log stand QMD (cm)",
    ba_metric   = "stand basal area (m2/ha)",
    bgi         = "stand mean bgi (climate/site driver); 2-piece spline knot at median",
    bgi_knot    = bgi_med,
    trt_decay   = "exp(-min(years_since_trt)/tau_m) over the stand",
    dstrb_decay = "exp(-min(years_since_dstrb)/tau_d) over the stand"),
  tau_m = TAU_M, tau_d = TAU_D,
  n_stand_records = nrow(ds),
  fixed_effects = fixed_list,
  sd_L1 = sd_L1,
  re_L1 = re_table,
  seed = 20260702,
  notes = "Stand mortality target M_stand feeds stand_constraint.py::stand_disaggregate_mortality to rescale summed tree hazards (proportional-hazard kappa)."
)
write_json(bundle, file.path(OUT_DIR, "stand_survival_bundle.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")

cat("\nStand-survival fixed effects (log-hazard scale):\n")
print(round(fx[, c("Estimate", "Est.Error", "Q2.5", "Q97.5"), drop = FALSE], 4))
cat("\nsd(1|L1):", round(sd_L1, 4), "\n")
gc(); cat("\nDONE_STAND_SURVIVAL\n")
