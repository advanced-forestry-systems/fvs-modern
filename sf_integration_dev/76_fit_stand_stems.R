#!/usr/bin/env Rscript
# =============================================================================
# 76_fit_stand_stems.R -- STAND-LEVEL García state-space STEM-DENSITY N(t)
# transition. There is NO fitted García-style stem-number transition in the
# repo (09_fit_stand_density.R fits only the Reineke SDIMAX self-thinning LINE;
# 35 is ingrowth/recruitment; 17 drives N purely by tree mortality). This script
# fits the minimal Bayesian state-space transition of stand stem number N so the
# tree-mortality constraint stand_constraint.py::stand_constrain_stems has a
# fitted N(t) TARGET:
#
#   N2 = surviving stems/ha at t2 (survivors of the N1 cohort; ingrowth excluded,
#        it belongs to the recruitment channel, not the survival channel this
#        constraint governs). deaths = N1 - N2 (>= 0).
#
#   cloglog exposure survival (SAME scale as 71_fit_stand_survival.R and the tree
#   mortality model), so the stem target is the count form of the survival target:
#     cloglog(M_stand) = log(H_stand) + log(YEARS) = eta_h + log(YEARS)
#     eta_h = b0 + b_lnN1*ln(N1) + b_topht*top_ht + b_rd*rd + b_lnqmd*ln_qmd + (1|EPA_L1)
#     H_stand = exp(eta_h);  S(T) = exp(-H_stand*T);  N2_target = N1 * S(T)
#
#   García state-space signature: next state N2 depends on the current state N1,
#   the driver TOP HEIGHT (self-thinning is height/size driven, not age), relative
#   density rd, and the interval. This is the density-transition analogue of the
#   García GADA top-height transition (75_fit_stand_topheight.R).
#
#   log(YEARS)-3.9 OFFSET TRICK: a cloglog exposure GLM with a raw log(YEARS)
#   offset can hit a log(0)/degenerate init failure; centering the exposure offset
#   at log(10) ~ 3.9 (median interval) keeps the linear predictor near 0 at init.
#   We use offset = log(YEARS) - 3.9 and add +3.9 back into the intercept at
#   evaluation time (folded into b0), so H_stand is unchanged.
#
# Usage:
#   Rscript 76_fit_stand_stems.R --n_sub=120000 --top_n=100 \
#     --pairs=<rds> --out_dir=<dir>
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
PAIRS <- ga("pairs", "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR <- ga("out_dir", "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_stems")
N_SUB <- as.integer(ga("n_sub", "120000"))
TOP_N_PER_HA <- as.numeric(ga("top_n", "100"))
OFFSET_CENTER <- as.numeric(ga("offset_center", "3.9"))   # log(YEARS) baseline offset (~log 10)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 76_fit_stand_stems.R ==\n top_n:", TOP_N_PER_HA, " n_sub:", N_SUB,
    " offset_center:", OFFSET_CENTER, "\n")

d <- trycatch_run(as.data.table(readRDS(PAIRS)), "read pairs")
if (is.null(d)) quit(status = 1)
cat(" pairs rows:", nrow(d), "\n")

# ---- row filter: live-at-t1 trees (same survival basis as 71) ----------------
if (!"TREESTATUS1" %in% names(d)) d[, TREESTATUS1 := 1L]
if (!"TPH_UNADJ1" %in% names(d) && "TPA1" %in% names(d)) d[, TPH_UNADJ1 := TPA1 * 2.4710538147]
d <- d[TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1, 2) &
       is.finite(TPH_UNADJ1) & TPH_UNADJ1 > 0 &
       is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
       is.finite(DBH1) & DBH1 >= 2.54]
cat(" live-at-t1 tree rows:", nrow(d), "\n")

d[, stand_key := paste(PLT_CN_cond1, CONDID_cond1, INVYR1, INVYR2, sep = "_")]
d[, died := as.numeric(TREESTATUS2 == 2)]

# top-height driver at t1 (tallest top_n_per_ha stems/ha) for the transition
top_ht <- function(ht, w, need) {
  ok <- is.finite(ht) & is.finite(w) & w > 0
  ht <- ht[ok]; w <- w[ok]
  if (!length(ht)) return(NA_real_)
  o <- order(-ht); w <- w[o]; ht <- ht[o]; cum <- cumsum(w)
  if (cum[length(cum)] <= need) return(sum(w * ht) / sum(w))
  full <- cum <= need; wt <- ifelse(full, w, 0)
  k <- which(!full)[1]; prev <- if (k > 1) cum[k - 1] else 0
  wt[k] <- need - prev; sum(wt * ht) / sum(wt)
}

agg <- d[, {
  N1 <- sum(TPH_UNADJ1)
  N2 <- sum(TPH_UNADJ1 * (1 - died))            # surviving stems/ha
  th <- if ("HT1" %in% names(.SD)) top_ht(HT1, TPH_UNADJ1, TOP_N_PER_HA) else NA_real_
  .(YEARS = YEARS[1], N1 = N1, N2 = N2, deaths = N1 - N2, top_ht = th,
    qmd = if ("QMD1" %in% names(.SD)) weighted.mean(QMD1, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    sdi_add = if ("sdi_additive1" %in% names(.SD)) weighted.mean(sdi_additive1, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    sdimax  = if ("SDImax_brms" %in% names(.SD)) weighted.mean(SDImax_brms, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    rd_row  = if ("rd_add" %in% names(.SD)) weighted.mean(rd_add, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    L1 = EPA_L1_CODE[1])
}, by = stand_key,
   .SDcols = intersect(c("HT1","QMD1","sdi_additive1","SDImax_brms","rd_add"), names(d))]
rm(d); gc()
cat(" stand-remeasurement stem records:", nrow(agg), "\n")

# ---- covariates + binomial response (deaths of N1) ---------------------------
agg[, rd := ifelse(is.finite(sdi_add) & is.finite(sdimax) & sdimax > 0, sdi_add / sdimax, rd_row)]
agg <- agg[is.finite(N1) & N1 > 0 & is.finite(N2) & N2 >= 0 & N2 <= N1 &
           is.finite(YEARS) & is.finite(rd) & rd > 0 & rd < 2]
agg[, ln_n1 := log(N1)]
agg[, ln_qmd := ifelse(is.finite(qmd) & qmd > 0, log(qmd), NA_real_)]
th_med <- suppressWarnings(median(agg$top_ht[is.finite(agg$top_ht)], na.rm = TRUE))
if (!is.finite(th_med)) th_med <- 15
agg[!is.finite(top_ht), top_ht := th_med]
agg[!is.finite(ln_qmd), ln_qmd := suppressWarnings(median(agg$ln_qmd[is.finite(agg$ln_qmd)], na.rm = TRUE))]
# log(YEARS) - 3.9 baseline-offset trick (avoid cloglog log(0) init failure)
agg[, log_years_off := log(YEARS) - OFFSET_CENTER]
agg[, L1 := as.character(L1)]
# integer binomial counts: trials = round(N1), deaths = round(N1 - N2)
agg[, trials_i := as.integer(round(pmax(N1, 1)))]
agg[, deaths_i := as.integer(pmin(pmax(round(N1 - N2), 0), trials_i))]
agg <- agg[!is.na(L1) & L1 != "" & is.finite(log_years_off) & is.finite(ln_n1) & is.finite(ln_qmd)]
cat(" modeling-ready records:", nrow(agg), "\n")
cat(" mean stand stem-mortality fraction:", round(mean(agg$deaths_i / agg$trials_i), 5), "\n")

# ---- subsample ---------------------------------------------------------------
n_use <- min(nrow(agg), N_SUB)
ds <- agg[sample(.N, n_use)]
cat(sprintf(" fit sample: %d of %d\n", nrow(ds), nrow(agg)))
rm(agg); gc()

# ---- Bayesian fit (brms): cloglog binomial + centered log(YEARS) offset -------
ok <- trycatch_run({ suppressPackageStartupMessages(library(brms)); TRUE }, "load brms")
if (is.null(ok)) quit(status = 1)

form <- bf(deaths_i | trials(trials_i) ~ ln_n1 + top_ht + rd + ln_qmd +
             (1 | L1) + offset(log_years_off))
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"))
fit <- trycatch_run(
  brm(form, data = ds, family = binomial(link = "cloglog"), prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, control = list(adapt_delta = 0.9),
      init = 0),   # init=0 + centered offset avoids log(0) init failure
  "brms stand-stems fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
# fold the OFFSET_CENTER back into the intercept so the reported eta_h is on the
# raw log(YEARS) exposure scale: H_stand = exp(eta_h); with centered offset the
# model intercept absorbed a -OFFSET_CENTER, so add it back for downstream use.
fx_adj <- fx
fx_adj["Intercept", "Estimate"] <- fx["Intercept", "Estimate"] - OFFSET_CENTER
sm <- as.data.frame(fx_adj); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, "stand_stems_fixed.csv"))

fe_names <- rownames(fx_adj)
fixed_list <- setNames(
  lapply(fe_names, function(nm) list(mean = unname(fx_adj[nm,"Estimate"]),
                                     sd = unname(fx[nm,"Est.Error"]),
                                     q2.5 = unname(fx[nm,"Q2.5"]) - if (nm=="Intercept") OFFSET_CENTER else 0,
                                     q97.5 = unname(fx[nm,"Q97.5"]) - if (nm=="Intercept") OFFSET_CENTER else 0)),
  fe_names)
re_table <- if (!is.null(re)) list(level = rownames(re), mean = unname(re[,"Estimate"]),
                                   sd = unname(re[,"Est.Error"])) else NULL
sd_L1 <- trycatch_run({ vc <- summary(fit)$random$L1; unname(vc["sd(Intercept)","Estimate"]) }, "sd L1")

bundle <- list(
  model = "stand-level García state-space STEM-DENSITY N(t) transition",
  scale = "cloglog binomial with centered log(YEARS)-3.9 exposure offset; intercept refolded to raw log(YEARS). Linear predictor eta_h = log stand hazard, SAME scale as 71_fit_stand_survival.R and the tree mortality model",
  response = "deaths_i | trials(trials_i); trials = round(N1 stems/ha), deaths = round(N1 - N2 survivors)",
  stand_key = "PLT_CN_cond1 x CONDID_cond1 x INVYR1 x INVYR2 (FIA plot-condition-remeasurement)",
  top_n_per_ha = TOP_N_PER_HA,
  offset_center = OFFSET_CENTER,
  transition = "H_stand = exp(eta_h); S(T)=exp(-H_stand*T); N2_target = N1 * S(T); eta_h = b0 + b_lnN1*ln(N1) + b_topht*top_ht + b_rd*rd + b_lnqmd*ln_qmd + z_L1",
  novelty = "No García-style stem-number transition existed in the repo (only the Reineke SDIMAX line in 09_fit_stand_density.R); this is the newly fitted N(t) state-space transition",
  covariates = list(
    ln_n1 = "log start stand stems/ha (current state N1)",
    top_ht = "stand top height (m), tallest top_n_per_ha stems/ha -- self-thinning driver",
    rd = "stand SDI / SDIMAX (sdi_additive1 / SDImax_brms)",
    ln_qmd = "log stand QMD (cm)"),
  n_stand_records = nrow(ds),
  fixed_effects = fixed_list,
  sd_L1 = sd_L1,
  re_L1 = re_table,
  seed = 20260702,
  notes = "N2 target feeds stand_constraint.py::stand_constrain_stems, which reuses the kappa proportional-hazard solve (N_target is the stems form of M_stand)."
)
write_json(bundle, file.path(OUT_DIR, "stand_stems_bundle.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")

cat("\nStand stem-density transition fixed effects (log-hazard scale, refolded):\n")
print(round(fx_adj[, c("Estimate","Est.Error","Q2.5","Q97.5"), drop = FALSE], 5))
cat("\nsd(1|L1):", round(sd_L1, 5), "\n")
gc(); cat("\nDONE_STAND_STEMS\n")
