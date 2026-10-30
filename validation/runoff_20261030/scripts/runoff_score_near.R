#!/usr/bin/env Rscript
## runoff_score_near.R OBS OUTDIR near_file1 [near_file2 ...]
## Near-term (single cycle = own remeasurement interval) accuracy, Paper 2 convention: bias, percent bias
## (100 x mean(pred - obs) / mean(obs)) and RMSE of BA (m2 ha-1), TPH (ha-1), QMD (cm) and common
## volume VOLC (m3 ha-1, 0.40 x BA x top height), stand-bootstrap 95% percentile CIs (2000 resamples,
## seed 20261028). Observed time 2 = survivor metrics from obs/nearterm_obs_frame.csv.
## Sets: "common" = stands that ran in EVERY arm of the variant (engine rows also seed_ok) = primary;
## "all" = every stand that ran in that arm (sensitivity). Also writes the year-0 seeding check.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); OBS <- a[1]; OUT <- a[2]; FILES <- a[-(1:2)]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NB <- 2000; set.seed(20261028)
ob <- fread(OBS)[obs_ok == TRUE]
P <- rbindlist(lapply(FILES, function(f) { d <- fread(f); d[, src := basename(f)]; d }), fill = TRUE)
P[, seed_ok := as.logical(seed_ok)]
P <- merge(P, ob[, .(plot_key, BA1o = BA1, TPH1o = TPH1, BAo = BA2, TPHo = TPH2, QMDo = QMD2, HTo = HTTOP2, VOLCo = VOLC2, BA_totalo = BA2_file, TPH_totalo = TPH2_file, k_obs = k_cycle, interval)], by = "plot_key")
cat("rows:", nrow(P), " rows whose run length differs from the observed interval (rounded):", P[years_run != k_obs, .N], "\n")
y0 <- P[, .(n = .N, seed_ok_share = mean(seed_ok %in% TRUE), med_ba0_over_balive = median(ba0_ft2ac / balive, na.rm = TRUE),
            med_BA0_over_obsBA1 = median(BA0 / BA1o, na.rm = TRUE), med_TPH0_over_obsTPH1 = median(TPH0 / TPH1o, na.rm = TRUE),
            med_k = as.numeric(median(k_obs))), by = .(variant, arm)][order(variant, arm)]
fwrite(y0, file.path(OUT, "near_year0_check.csv")); print(y0)
P[, ok := is.finite(BA) & is.finite(TPH) & TPH > 0 & (seed_ok %in% TRUE)]
narm <- P[, .(narm = uniqueN(arm)), by = variant]
cs <- merge(P[ok == TRUE, .(na = uniqueN(arm)), by = .(variant, traj)], narm, by = "variant")[na == narm, .(variant, traj)]
P[, common := FALSE]; P[cs, common := TRUE, on = .(variant, traj)]
cat("common stand set by variant:\n"); print(P[common == TRUE, .(stands = uniqueN(traj)), by = variant])
P[, `:=`(BA_total = BA, TPH_total = TPH)]   # sensitivity: against the file BA2/TPH2 (all live at time 2, ingrowth included)
mets <- c("BA", "TPH", "QMD", "VOLC", "BA_total", "TPH_total")
st <- function(p, o) c(bias = mean(p - o), pbias = 100 * mean(p - o) / mean(o), rmse = sqrt(mean((p - o)^2)))
res <- list()
for (set in c("common", "all")) {
  Q <- if (set == "common") P[common == TRUE] else P[ok == TRUE]
  for (vv in sort(unique(Q$variant))) for (ar in sort(unique(Q[variant == vv]$arm))) {
    d <- Q[variant == vv & arm == ar]
    for (m in mets) { p <- d[[m]]; o <- d[[paste0(m, "o")]]; k <- is.finite(p) & is.finite(o); p <- p[k]; o <- o[k]; if (length(p) < 2) next
      est <- st(p, o); bs <- replicate(NB, { i <- sample.int(length(p), replace = TRUE); st(p[i], o[i]) })
      ci <- apply(bs, 1, quantile, c(.025, .975), na.rm = TRUE)
      res[[length(res) + 1]] <- data.table(set = set, variant = vv, arm = ar, metric = m, n = length(p), obs_mean = mean(o), pred_mean = mean(p),
        bias = est[["bias"]], bias_lo = ci[1, "bias"], bias_hi = ci[2, "bias"], pbias = est[["pbias"]], pbias_lo = ci[1, "pbias"], pbias_hi = ci[2, "pbias"],
        rmse = est[["rmse"]], rmse_lo = ci[1, "rmse"], rmse_hi = ci[2, "rmse"]) } } }
R <- rbindlist(res); fwrite(R, file.path(OUT, "near_scores.csv"))
print(R[set == "common", .(variant, arm, metric, n, obs = round(obs_mean, 1), pred = round(pred_mean, 1), bias = round(bias, 2),
  pbias = round(pbias, 1), lo = round(pbias_lo, 1), hi = round(pbias_hi, 1), rmse = round(rmse, 2))], nrows = 500)
fwrite(P, file.path(OUT, "near_joined.csv"))
cat("NEAR_SCORE_DONE\n")
