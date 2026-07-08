# 62k: per-species best-AIC HT-DBH model selection -> single consistent coefficient table
# Aaron decision 2026-07-07: per-species lowest-AIC among acceptable fits, all 6 forms dispatched by model_id.
suppressWarnings(suppressMessages(library(data.table)))
dir <- "/fs/scratch/PUOM0008/crsfaaron/wt-htdbh/sf_integration_dev/htdbh_native/marshall_fits"
mods <- rbindlist(lapply(1:6, function(m) fread(file.path(dir, sprintf("CONUS_HDfit_Model_%d.CSV", m)))), fill=TRUE)
setnames(mods, tolower(names(mods)))
# acceptable tiers: 0 = clean converged, right sign, significant; 4 = converged/right-sign but some nonsig param
mods[, tier := fifelse(ierr==0, 1L, fifelse(ierr==4, 2L, NA_integer_))]
acc <- mods[!is.na(tier)]
setorder(acc, spcd, tier, aic)                 # prefer tier 1, then lowest AIC
best <- acc[, .SD[1], by=spcd]
best[, `:=`(model_id=model, B1=b1, B2=b2, B3=b3)]
out <- best[, .(SPCD=spcd, model_id, B1, B2, B3, nobs, AIC=aic, RMSE=rmse, meanBIAS=meanbias, ierr, tier)]
setorder(out, SPCD)
fwrite(out, "/fs/scratch/PUOM0008/crsfaaron/wt-htdbh/sf_integration_dev/htdbh_native/conus_htdbh_coefficients.csv")
allspp <- fread(file.path(dir,"CONUS_HDfit_Data.CSV")); setnames(allspp, tolower(names(allspp)))
nofit <- setdiff(allspp$spcd, out$SPCD)
cat("=== total species with data:", nrow(allspp), " | selected:", nrow(out), " | native-fallback (no acceptable fit):", length(nofit), "===\n")
cat("=== model-mix (per-species best-AIC) ===\n"); print(out[, .N, by=model_id][order(model_id)])
cat("=== tier mix (1=clean,2=nonsig-param accepted) ===\n"); print(out[, .N, by=tier][order(tier)])
cat("=== native-fallback SPCDs ===\n"); print(sort(nofit))
cat("=== head of coefficient table ===\n"); print(head(out, 6))
