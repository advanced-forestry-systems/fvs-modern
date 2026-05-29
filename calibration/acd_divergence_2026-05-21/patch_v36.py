#!/usr/bin/env python3
# Build cardinal_acadgy_v33v34test_v36.R = v32 harness scaled up to 300 plots
# with seed=2028 (different from v30 seed=42 and v32 seed=2027). Sources v33
# coefficients via apply_density_correction.R, runs 12.3.9 production posture,
# applies stand-level + tree-level reconciliation, reports diagnostics. Larger
# fresh sample = true out-of-sample test of generalization.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_v31test_v32.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_v33v34test_v36.R"
s = open(src).read()

s = s.replace("[v32]", "[v36]")
s = s.replace("[v30]", "[v36]")
s = s.replace("acadgy_v31test_v32_results.csv", "acadgy_v33v34test_v36_results.csv")
s = s.replace("acdgy_v31test_v32_perplot.csv", "acdgy_v33v34test_v36_perplot.csv")

# Up N_PLOTS to 300 and seed to 2028
s = s.replace("N_PLOTS <- 100L", "N_PLOTS <- 300L")
s = s.replace("set.seed(2027)", "set.seed(2028)")

# The v32 script had a crash on the multivariate regression. Replace the
# entire tail (residual analysis + v31 application) with a clean v33+v34 test.
import re
pat = re.compile(r'# ---- v31 out-of-sample test ----.*?cat\("\\n=== Done\. Fresh-sample perplot at acdgy_v31test_v32_perplot\.csv ===\\n"\)',
                 re.DOTALL)
new_tail = (
    '# ---- v33+v34 out-of-sample test on 300-plot fresh sample ----\n'
    'source("/users/PUOM0008/crsfaaron/acadgy_fia_verify/apply_density_correction.R")\n'
    'cat(sprintf("[v36] v33 coefficients: a=%.4f b=%.6f lower=%.0f upper=%.0f (fit on n=%d)\\n",\n'
    '    ACD_DENSITY_CORRECTION$a, ACD_DENSITY_CORRECTION$b,\n'
    '    ACD_DENSITY_CORRECTION$lower_cap, ACD_DENSITY_CORRECTION$upper_cap,\n'
    '    ACD_DENSITY_CORRECTION$n))\n'
    '\n'
    'pp <- get("perplot_csi_scale_0.7")\n'
    'pp$PLT_CN <- as.character(samp$PLT_CN_t1_str)[seq_len(nrow(pp))]\n'
    'pp$BA_t1 <- samp$BA_t1[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp$interval <- samp$interval_years[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp_complete <- pp[complete.cases(pp[, c("BA_pred","BA_t2","BA_t1","TPA","TPA_t2")]), ]\n'
    'cat(sprintf("[v36] n complete cases = %d / 300\\n", nrow(pp_complete)))\n'
    '\n'
    '# Stand-level v33\n'
    'pp_complete$BA_corrected <- apply_density_correction(pp_complete$BA_pred, pp_complete$BA_t1)\n'
    '\n'
    '# Tree-level v34 simulated: scale TPA by the same scale factor (since we do not\n'
    '# have tree-level data in pp, we simulate the effect by scaling TPA in proportion)\n'
    'pp_complete$BA_t1_safe <- pmin(pmax(pp_complete$BA_t1, 0), 400)\n'
    'raw_corr  <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * pp_complete$BA_t1_safe\n'
    'bnd_corr  <- pmax(ACD_DENSITY_CORRECTION$lower_cap, pmin(ACD_DENSITY_CORRECTION$upper_cap, raw_corr))\n'
    'raw_scale <- ifelse(pp_complete$BA_pred > 0, (pp_complete$BA_pred - bnd_corr) / pp_complete$BA_pred, 1)\n'
    'bnd_scale <- pmax(0.7, pmin(1.0, raw_scale))\n'
    'pp_complete$BA_treerecon <- pp_complete$BA_pred * bnd_scale\n'
    'pp_complete$TPA_treerecon <- pp_complete$TPA * bnd_scale\n'
    '\n'
    'bias <- function(y, yhat) 100*(mean(yhat) - mean(y))/mean(y)\n'
    'r2   <- function(y, yhat) 1 - sum((y-yhat)^2)/sum((y-mean(y))^2)\n'
    '\n'
    'cat("\\n=== v36 out-of-sample test (300-plot fresh ME FIA sample, seed=2028) ===\\n")\n'
    'cat(sprintf("Uncorrected:           BA bias = %+.2f%%   R^2 = %.4f   TPA bias = %+.2f%%\\n",\n'
    '    bias(pp_complete$BA_t2, pp_complete$BA_pred), r2(pp_complete$BA_t2, pp_complete$BA_pred),\n'
    '    bias(pp_complete$TPA_t2, pp_complete$TPA)))\n'
    'cat(sprintf("v33 stand-level corr:  BA bias = %+.2f%%   R^2 = %.4f   TPA bias = %+.2f%% (TPA unchanged)\\n",\n'
    '    bias(pp_complete$BA_t2, pp_complete$BA_corrected), r2(pp_complete$BA_t2, pp_complete$BA_corrected),\n'
    '    bias(pp_complete$TPA_t2, pp_complete$TPA)))\n'
    'cat(sprintf("v34 tree-level recon:  BA bias = %+.2f%%   R^2 = %.4f   TPA bias = %+.2f%% (TPA scaled)\\n",\n'
    '    bias(pp_complete$BA_t2, pp_complete$BA_treerecon), r2(pp_complete$BA_t2, pp_complete$BA_treerecon),\n'
    '    bias(pp_complete$TPA_t2, pp_complete$TPA_treerecon)))\n'
    '\n'
    'cat("\\n=== Per BA_t1 quartile diagnostic ===\\n")\n'
    'pp_complete$q <- cut(pp_complete$BA_t1,\n'
    '   breaks=quantile(pp_complete$BA_t1, probs=seq(0,1,0.25), na.rm=TRUE),\n'
    '   include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4"))\n'
    'qtbl <- do.call(rbind, lapply(split(pp_complete, pp_complete$q), function(g) data.frame(\n'
    '   n=nrow(g), BA_t1_mean=round(mean(g$BA_t1, na.rm=TRUE), 1),\n'
    '   raw_BA_bias=round(bias(g$BA_t2, g$BA_pred), 2),\n'
    '   v33_BA_bias=round(bias(g$BA_t2, g$BA_corrected), 2),\n'
    '   v34_BA_bias=round(bias(g$BA_t2, g$BA_treerecon), 2),\n'
    '   v34_TPA_bias=round(bias(g$TPA_t2, g$TPA_treerecon), 2),\n'
    '   mean_scale=round(mean(bnd_scale[pp_complete$q == names(table(pp_complete$q))[1]]), 3))))\n'
    'qtbl$q <- rownames(qtbl)\n'
    'print(format(qtbl, digits=3))\n'
    '\n'
    'write.csv(pp_complete, file.path(OUT_DIR, "acdgy_v33v34test_v36_perplot.csv"), row.names=FALSE)\n'
    'cat(sprintf("\\n=== Done. Per-plot CSV at acdgy_v33v34test_v36_perplot.csv (n=%d) ===\\n", nrow(pp_complete)))\n'
)
m = pat.search(s)
if m:
    s = pat.sub(new_tail, s, count=1)
else:
    s = s + "\n" + new_tail + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "n_plots_ok", "N_PLOTS <- 300L" in s,
      "seed_ok", "set.seed(2028)" in s,
      "v33_apply_ok", "apply_density_correction(" in s,
      "v34_recon_ok", "BA_treerecon" in s)
