#!/usr/bin/env python3
# Build cardinal_acadgy_residualcal_v30.R = v28 (12.3.9 CSI_SCALE=0.7 harness)
# but emit per-plot predictions, not just aggregates. Then fit a residual
# regression locally: residual = f(BGI, SICOND, FORTYPCD, BAL_t1) to see if any
# carry signal that a post-hoc calibration could exploit. The v27c/v29 input-
# based tests were null; this asks the orthogonal question via the residual.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_csiscale_v28.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_residualcal_v30.R"
s = open(src).read()

s = s.replace("[v28]", "[v30]")
s = s.replace("acadgy_csiscale_v28_results.csv", "acadgy_residualcal_v30_results.csv")

# Modify summ() to also retain the per-plot dataframe in a global list
old_summ = ("summ <- function(tag, d) {\n"
            "  o<-d$BA_t2; p<-d$BA_pred; r2<-1-sum((p-o)^2)/sum((o-mean(o))^2)\n"
            "  cat(sprintf(\"%-14s BA_bias=%+.1f%% R2=%.3f TPA=%.0f(obs %.0f) QMD=%.2f(obs %.2f)\\n\",\n"
            "              tag, 100*(mean(p)-mean(o))/mean(o), r2, mean(d$TPA,na.rm=T), mean(d$TPA_t2,na.rm=T), mean(d$QMD_in,na.rm=T), mean(d$QMD_t2,na.rm=T)))\n"
            "  data.frame(config=tag, BA_obs=mean(o), BA_pred=mean(p), BA_bias_pct=100*(mean(p)-mean(o))/mean(o), BA_r2=r2,\n"
            "             TPA=mean(d$TPA,na.rm=T), TPA_obs=mean(d$TPA_t2,na.rm=T), QMD=mean(d$QMD_in,na.rm=T), QMD_obs=mean(d$QMD_t2,na.rm=T))\n"
            "}")
new_summ = ("summ <- function(tag, d) {\n"
            "  o<-d$BA_t2; p<-d$BA_pred; r2<-1-sum((p-o)^2)/sum((o-mean(o))^2)\n"
            "  cat(sprintf(\"%-14s BA_bias=%+.1f%% R2=%.3f TPA=%.0f(obs %.0f) QMD=%.2f(obs %.2f)\\n\",\n"
            "              tag, 100*(mean(p)-mean(o))/mean(o), r2, mean(d$TPA,na.rm=T), mean(d$TPA_t2,na.rm=T), mean(d$QMD_in,na.rm=T), mean(d$QMD_t2,na.rm=T)))\n"
            "  d$config <- tag\n"
            "  d$PLT_CN <- as.character(samp$PLT_CN_t1_str)[seq_len(nrow(d))]\n"
            "  assign(paste0(\"perplot_\", tag), d, envir=.GlobalEnv)\n"
            "  data.frame(config=tag, BA_obs=mean(o), BA_pred=mean(p), BA_bias_pct=100*(mean(p)-mean(o))/mean(o), BA_r2=r2,\n"
            "             TPA=mean(d$TPA,na.rm=T), TPA_obs=mean(d$TPA_t2,na.rm=T), QMD=mean(d$QMD_in,na.rm=T), QMD_obs=mean(d$QMD_t2,na.rm=T))\n"
            "}")
assert old_summ in s, "summ function not found"
s = s.replace(old_summ, new_summ, 1)

# Replace run block with two configs (baseline + production), then residual
# analysis with BGI, SICOND, FORTYPCD, BAL_t1
import re
pat = re.compile(r'rows <- list\(\)\s*\ncat\("\[v30\] csi_scale_1\.0.*?print\(format\(res, digits=4\)\)',
                 re.DOTALL)
new_block = (
    'rows <- list()\n'
    'cat("[v30] csi_scale_0.7  (12.3.9, MORTCAL on, CutPoint 0, CSI_SCALE = 0.7; production posture)\\n")\n'
    'rows[["a"]] <- summ("csi_scale_0.7",  run_cfg(TRUE, "Y", 0, 0.7))\n'
    'cat("[v30] csi_scale_1.0  (12.3.9, MORTCAL on, CutPoint 0, no CSI_SCALE; baseline)\\n")\n'
    'rows[["b"]] <- summ("csi_scale_1.0",  run_cfg(TRUE, "Y", 0, NULL))\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_residualcal_v30_results.csv"), row.names=FALSE)\n'
    'print(format(res, digits=4))\n'
    '\n'
    '# ---- Per-plot residual analysis ----\n'
    '# Join per-plot predictions with BGI and vdat covariates, regress residual.\n'
    '# Residual is in BA_pred - BA_t2 (positive = model overshoots, negative =\n'
    '# undershoots). If BGI explains residual variance, a post-hoc calibration\n'
    '# factor is possible.\n'
    'pp <- get("perplot_csi_scale_0.7")\n'
    'pp$residual <- pp$BA_pred - pp$BA_t2\n'
    'pp$relresid <- pp$residual / pp$BA_t2\n'
    'bgi_csv <- read.csv("/users/PUOM0008/crsfaaron/acadgy_fia_verify/me_bgi_by_pltcn.csv",\n'
    '                    colClasses=c("CN"="character"))\n'
    'bgi_csv$BGI[is.na(bgi_csv$BGI) | bgi_csv$BGI <= 0] <- NA\n'
    'pp$BGI <- setNames(bgi_csv$BGI, bgi_csv$CN)[as.character(pp$PLT_CN)]\n'
    'pp$SICOND <- samp$SICOND[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp$FORTYPCD <- samp$FORTYPCD[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp$BA_t1 <- samp$BA_t1[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp$interval <- samp$interval_years[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp$ClimateSI_ft <- samp$ClimateSI_ft[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp_complete <- pp[complete.cases(pp[, c("residual","BGI","SICOND","BA_t1","ClimateSI_ft")]), ]\n'
    'write.csv(pp_complete, file.path(OUT_DIR, "acadgy_residualcal_v30_perplot.csv"), row.names=FALSE)\n'
    'cat(sprintf("\\n=== v30 per-plot residual analysis: n=%d complete cases ===\\n", nrow(pp_complete)))\n'
    '\n'
    '# Univariate regressions on residual\n'
    'forms <- list(\n'
    '  "residual ~ BGI"               = residual ~ BGI,\n'
    '  "residual ~ log(BGI)"          = residual ~ log(BGI),\n'
    '  "residual ~ I(BGI^2)"          = residual ~ I(BGI^2),\n'
    '  "relresid ~ BGI"               = relresid ~ BGI,\n'
    '  "relresid ~ log(BGI)"          = relresid ~ log(BGI),\n'
    '  "residual ~ SICOND"            = residual ~ SICOND,\n'
    '  "residual ~ ClimateSI_ft"      = residual ~ ClimateSI_ft,\n'
    '  "residual ~ BA_t1"             = residual ~ BA_t1,\n'
    '  "residual ~ interval"          = residual ~ interval,\n'
    '  "residual ~ BGI + SICOND + BA_t1" = residual ~ BGI + SICOND + BA_t1,\n'
    '  "residual ~ BGI * BA_t1"       = residual ~ BGI * BA_t1\n'
    ')\n'
    'for (label in names(forms)) {\n'
    '  m <- lm(forms[[label]], data=pp_complete)\n'
    '  cat(sprintf("[v30] %-40s R^2 = %.4f  p(F) = %.4g\\n",\n'
    '              label, summary(m)$r.squared,\n'
    '              pf(summary(m)$fstatistic[1], summary(m)$fstatistic[2], summary(m)$fstatistic[3], lower.tail=FALSE)))\n'
    '}\n'
    '\n'
    '# Forest type marginal effects on residual (mean residual by FORTYPCD bucket)\n'
    'pp_complete$FT_bucket <- cut(pp_complete$FORTYPCD,\n'
    '   breaks=c(0,400,500,800,1000), labels=c("Spruce-fir","Mixedwood","Hardwood","Other"))\n'
    'cat("\\n=== Mean residual by forest type ===\\n")\n'
    'print(aggregate(residual ~ FT_bucket, data=pp_complete, FUN=function(x) c(mean=mean(x), n=length(x))))\n'
    '\n'
    'cat("\\n=== Done. Per-plot CSV at acadgy_residualcal_v30_perplot.csv ===\\n")'
)
m = pat.search(s)
if m:
    s = pat.sub(new_block, s, count=1)
else:
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "perplot_ok", "perplot_csi_scale_0.7" in s,
      "regression_ok", "residual ~ BGI" in s,
      "fortype_ok", "FT_bucket" in s)
