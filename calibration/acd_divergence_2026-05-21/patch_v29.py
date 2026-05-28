#!/usr/bin/env python3
# Build cardinal_acadgy_bgicsi_v29.R = v28 (12.3.9 CSI_SCALE harness) with two
# additional configs that route the ME BGI raster through ops$CSI_SCALE per
# stand. v27c showed BGI as a per-tree dDBH.mult is null; this test asks
# whether routing the same BGI through the climate channel (where CSI feeds
# dDBH AND ingrowth AND height) reveals plot-level signal that the per-tree
# version missed.
#
# Configs:
#   csi_scale_1.0      v24 sanity (CSI_SCALE not set)
#   csi_scale_0.7      v25 production recommendation (uniform 0.7)
#   csi_bgi_recenter   ops$CSI_SCALE = BGI / mean(BGI) per stand
#   csi_bgi_07x        ops$CSI_SCALE = 0.7 * (BGI / mean(BGI)) per stand
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_csiscale_v28.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_bgicsi_v29.R"
s = open(src).read()

s = s.replace("[v28]", "[v29]")
s = s.replace("acadgy_csiscale_v28_results.csv", "acadgy_bgicsi_v29_results.csv")

# Inject BGI load right after stand_init is built
anchor = "ops0 <- list(verbose = FALSE, INGROWTH = \"Y\", MinDBH = 3.0)"
assert anchor in s
inject = (
    "# v29: BGI per-stand CSI_SCALE. Load BGI extracted via gdallocationinfo from\n"
    "# ME_BGI_V1.tif. Recenter so mean of matched plots equals 1.0 to preserve\n"
    "# mean climate weighting and let BGI variance modulate per-stand CSI.\n"
    "bgi_csv <- read.csv(\"/users/PUOM0008/crsfaaron/acadgy_fia_verify/me_bgi_by_pltcn.csv\",\n"
    "                    colClasses=c(\"CN\"=\"character\"))\n"
    "bgi_csv$BGI[is.na(bgi_csv$BGI) | bgi_csv$BGI <= 0] <- NA\n"
    "bgi_lookup <- setNames(bgi_csv$BGI, bgi_csv$CN)\n"
    "bgi_per_plot <- bgi_lookup[as.character(samp$PLT_CN_t1_str)]\n"
    "bgi_mean_pos <- mean(bgi_per_plot, na.rm=TRUE)\n"
    "bgi_per_plot[is.na(bgi_per_plot)] <- bgi_mean_pos\n"
    "bgi_factor   <- setNames(bgi_per_plot / bgi_mean_pos, samp$PLT_CN_t1_str)\n"
    "cat(sprintf(\"[v29] BGI loaded: n=%d, mean=%.1f, factor range=[%.3f,%.3f]\\n\",\n"
    "    sum(!is.na(bgi_lookup[as.character(samp$PLT_CN_t1_str)])),\n"
    "    bgi_mean_pos, min(bgi_factor), max(bgi_factor)))\n"
    "\n"
    + anchor
)
s = s.replace(anchor, inject, 1)

# Extend p1y to also accept a per-stand csi_scale lookup
old_p1y = ("p1y <- function(trees, mortcal, ingrowth=\"Y\", cutpoint=NULL, csi_scale=NULL) {\n"
           "  pc <- list()\n"
           "  for (sid in unique(trees$STAND)) {\n"
           "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
           "    if (nrow(sub) == 0) next\n"
           "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
           "    if (!is.null(csi_scale)) ops$CSI_SCALE <- csi_scale\n"
           "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }")
new_p1y = ("p1y <- function(trees, mortcal, ingrowth=\"Y\", cutpoint=NULL, csi_scale=NULL, csi_per_stand=NULL) {\n"
           "  pc <- list()\n"
           "  for (sid in unique(trees$STAND)) {\n"
           "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
           "    if (nrow(sub) == 0) next\n"
           "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
           "    if (!is.null(csi_scale)) ops$CSI_SCALE <- csi_scale\n"
           "    if (!is.null(csi_per_stand)) {\n"
           "      f <- csi_per_stand[as.character(sid)]\n"
           "      if (!is.na(f) && is.finite(f)) {\n"
           "        if (!is.null(ops$CSI_SCALE)) ops$CSI_SCALE <- ops$CSI_SCALE * as.numeric(f)\n"
           "        else                          ops$CSI_SCALE <- as.numeric(f)\n"
           "      }\n"
           "    }\n"
           "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }")
assert old_p1y in s
s = s.replace(old_p1y, new_p1y, 1)

# Update run_cfg signature
s = s.replace(
    'run_cfg <- function(mortcal, ingrowth="Y", cutpoint=NULL, csi_scale=NULL) {',
    'run_cfg <- function(mortcal, ingrowth="Y", cutpoint=NULL, csi_scale=NULL, csi_per_stand=NULL) {', 1)
s = s.replace(
    'cur <- p1y(cur, mortcal, ingrowth, cutpoint, csi_scale); if (is.null(cur)) break',
    'cur <- p1y(cur, mortcal, ingrowth, cutpoint, csi_scale, csi_per_stand); if (is.null(cur)) break', 1)

# Replace the run block with four configs
import re
pat = re.compile(r'rows <- list\(\)\s*\ncat\("\[v29\] csi_scale_1\.0.*?print\(format\(res, digits=4\)\)',
                 re.DOTALL)
new_block = (
    'rows <- list()\n'
    'cat("[v29] csi_scale_1.0   (12.3.9, no CSI_SCALE; sanity match v24)\\n")\n'
    'rows[["a"]] <- summ("csi_scale_1.0",   run_cfg(TRUE, "Y", 0, NULL, NULL))\n'
    'cat("[v29] csi_scale_0.7   (12.3.9, CSI_SCALE = 0.7 uniform)\\n")\n'
    'rows[["b"]] <- summ("csi_scale_0.7",   run_cfg(TRUE, "Y", 0, 0.7, NULL))\n'
    'cat("[v29] csi_bgi_recenter(12.3.9, CSI_SCALE = BGI / mean(BGI) per stand; preserves mean climate)\\n")\n'
    'rows[["c"]] <- summ("csi_bgi_recenter",run_cfg(TRUE, "Y", 0, NULL, bgi_factor))\n'
    'cat("[v29] csi_bgi_0.7x    (12.3.9, CSI_SCALE = 0.7 * BGI / mean(BGI); stacked)\\n")\n'
    'rows[["d"]] <- summ("csi_bgi_0.7x",    run_cfg(TRUE, "Y", 0, 0.7, bgi_factor))\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_bgicsi_v29_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA: does BGI variance via CSI_SCALE help where per-tree dDBH.mult did not? ===\\n")\n'
    'print(format(res, digits=4))'
)
m = pat.search(s)
if m:
    s = pat.sub(new_block, s, count=1)
else:
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "csi_per_stand_ok", "csi_per_stand" in s,
      "bgi_factor_ok", "bgi_factor" in s,
      "configs_ok", "csi_bgi_recenter" in s and "csi_bgi_0.7x" in s)
