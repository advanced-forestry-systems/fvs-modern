#!/usr/bin/env python3
# Build cardinal_acadgy_v31test_v32.R = v30 harness but seeded differently to
# draw a fresh 100-plot ME FIA sample (no overlap guaranteed but different
# selection). Sources v31's apply_density_correction.R, runs production posture
# (12.3.9 + CSI_SCALE=0.7), applies v31 coefficients WITHOUT refitting, and
# reports test-set bias and R^2. If the result matches the v31 holdout
# (~+2 percent bias, ~0.48 corrected R^2), v31 generalizes. If it degrades
# significantly, the v30 fit overcaught noise.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_residualcal_v30.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_v31test_v32.R"
s = open(src).read()

# Tag swap
s = s.replace("[v30]", "[v32]")
s = s.replace("acadgy_residualcal_v30_results.csv", "acadgy_v31test_v32_results.csv")
s = s.replace("acadgy_residualcal_v30_perplot.csv", "acdgy_v31test_v32_perplot.csv")

# Change seed for fresh sample
s = s.replace("set.seed(42)", "set.seed(2027)")

# Replace the v30 regression block with v31 application + out-of-sample stats
import re
# Find the regression analysis block (from "# Per-plot residual analysis" through end of file or print(format(res...))
pat = re.compile(
    r'# ---- Per-plot residual analysis ----.*?cat\("\\n=== Done\. Per-plot CSV at acadgy_residualcal_v30_perplot\.csv ===\\n"\)',
    re.DOTALL)
new_tail = (
    '# ---- v31 out-of-sample test ----\n'
    '# Source v31 correction coefficients and apply WITHOUT refit. Then compare\n'
    '# uncorrected vs corrected bias and R^2 on this fresh sample.\n'
    'source("/users/PUOM0008/crsfaaron/acadgy_fia_verify/apply_density_correction.R")\n'
    'cat(sprintf("v31 coefficients: a=%.4f b=%.6f cap=%.0f (fit on n=%d, holdout R^2 mean=%.3f)\\n",\n'
    '    ACD_DENSITY_CORRECTION$a, ACD_DENSITY_CORRECTION$b,\n'
    '    ACD_DENSITY_CORRECTION$cap, ACD_DENSITY_CORRECTION$n,\n'
    '    ACD_DENSITY_CORRECTION$r2_holdout_mean))\n'
    '\n'
    'pp <- get("perplot_csi_scale_0.7")\n'
    'pp$PLT_CN <- as.character(samp$PLT_CN_t1_str)[seq_len(nrow(pp))]\n'
    'pp$BA_t1 <- samp$BA_t1[match(pp$PLT_CN, samp$PLT_CN_t1_str)]\n'
    'pp_complete <- pp[complete.cases(pp[, c("BA_pred","BA_t2","BA_t1")]), ]\n'
    'cat(sprintf("[v32] n complete cases = %d\\n", nrow(pp_complete)))\n'
    '\n'
    'pp_complete$BA_corrected <- apply_density_correction(pp_complete$BA_pred, pp_complete$BA_t1)\n'
    '\n'
    'bias <- function(y, yhat) 100*(mean(yhat) - mean(y))/mean(y)\n'
    'r2   <- function(y, yhat) 1 - sum((y-yhat)^2)/sum((y-mean(y))^2)\n'
    '\n'
    'raw_bias    <- bias(pp_complete$BA_t2, pp_complete$BA_pred)\n'
    'raw_r2      <- r2(pp_complete$BA_t2, pp_complete$BA_pred)\n'
    'corr_bias   <- bias(pp_complete$BA_t2, pp_complete$BA_corrected)\n'
    'corr_r2     <- r2(pp_complete$BA_t2, pp_complete$BA_corrected)\n'
    '\n'
    'cat("\\n=== v32 fresh-sample out-of-sample validation of v31 ===\\n")\n'
    'cat(sprintf("Uncorrected: BA_bias = %+.2f%%   R^2 = %.4f\\n", raw_bias, raw_r2))\n'
    'cat(sprintf("v31 corrected: BA_bias = %+.2f%%   R^2 = %.4f\\n", corr_bias, corr_r2))\n'
    'cat(sprintf("Closure: %.2f pp of bias  |  %+.4f R^2 lift\\n", raw_bias - corr_bias, corr_r2 - raw_r2))\n'
    'cat(sprintf("\\nExpected (v31 holdout): bias closes from +10.86%% to +2.11%%; R^2 lifts 0.379 -> 0.479\\n"))\n'
    '\n'
    'write.csv(pp_complete, file.path(OUT_DIR, "acadgy_v31test_v32_perplot.csv"), row.names=FALSE)\n'
    '\n'
    '# Also report by BA_t1 quartile to see if the correction is uniform across\n'
    '# density regimes (a key generalization risk check).\n'
    'pp_complete$BA_t1_q <- cut(pp_complete$BA_t1, breaks=quantile(pp_complete$BA_t1, probs=seq(0,1,0.25), na.rm=TRUE),\n'
    '                            include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4"))\n'
    'cat("\\n=== Per-quartile diagnostic ===\\n")\n'
    'q_table <- do.call(rbind, lapply(split(pp_complete, pp_complete$BA_t1_q), function(g) data.frame(\n'
    '  n=nrow(g), BA_t1_mean=mean(g$BA_t1, na.rm=TRUE),\n'
    '  raw_bias=bias(g$BA_t2, g$BA_pred), corr_bias=bias(g$BA_t2, g$BA_corrected))))\n'
    'q_table <- cbind(quartile=rownames(q_table), q_table)\n'
    'print(format(q_table, digits=3))\n'
    '\n'
    'cat("\\n=== Done. Fresh-sample perplot at acadgy_v31test_v32_perplot.csv ===\\n")'
)
m = pat.search(s)
if m:
    s = pat.sub(new_tail, s, count=1)
else:
    s = s + "\n" + new_tail + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "seed_ok", "set.seed(2027)" in s,
      "v31_apply_ok", "apply_density_correction(" in s,
      "quartile_ok", "Per-quartile diagnostic" in s)
