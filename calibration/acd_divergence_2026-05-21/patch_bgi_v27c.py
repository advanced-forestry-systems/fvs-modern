#!/usr/bin/env python3
# Build cardinal_acadgy_bgi_v27c.R = v27b + two more configs that use
# ME_BGI_V1.tif extracted at each plot's lat/lon as a plot-level dDBH.mult.
# Extraction: project FIA WGS84 lat/lon (EPSG:4326) to NAD83 UTM 19N
# (EPSG:26919), then terra::extract.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_siteindex_v27b.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_bgi_v27c.R"
s = open(src).read()

s = s.replace("[v27b]", "[v27c]")
s = s.replace("acadgy_siteindex_v27b_results.csv", "acadgy_bgi_v27c_results.csv")
s = s.replace("siteindex_v27b", "bgi_v27c")

# Inject BGI raster extraction right after stand_init_fvssi is built. Use terra
# (lighter than raster, available in R 4.4 stack). Match plots by PLT_CN_t1.
anchor = ("stand_init_fvssi <- unique(stand_init_fvssi)\n"
          "cat(sprintf(\"[v27c] CSI mean=%.2f m | SICOND mean=%.2f m | FVS_SITE_INDEX mean=%.2f m\\n\",\n"
          "    mean(stand_init_csi$CSI, na.rm=TRUE),\n"
          "    mean(stand_init_sicond$CSI, na.rm=TRUE),\n"
          "    mean(stand_init_fvssi$CSI, na.rm=TRUE)))")
assert anchor in s, "v27b stand_init anchor not found"
extract_block = anchor + (
    "\n\n"
    "# v27c: extract ME BGI raster at each plot's lat/lon, store per-stand.\n"
    "# ME_PLOT.csv has LAT, LON in WGS84 decimal degrees. Raster is NAD83 UTM 19N.\n"
    "if (requireNamespace(\"terra\", quietly=TRUE)) {\n"
    "  bgi_rast <- terra::rast(\"/users/PUOM0008/crsfaaron/raster_layers/bgi/ME_BGI_V1.tif\")\n"
    "  bgi_proj <- terra::crs(bgi_rast)\n"
    "  me_plot_xy <- me_plot[!is.na(me_plot$LAT) & !is.na(me_plot$LON) & me_plot$CN %in% samp$PLT_CN_t1_str, ]\n"
    "  pts_ll <- terra::vect(cbind(me_plot_xy$LON, me_plot_xy$LAT), type=\"points\", crs=\"EPSG:4326\")\n"
    "  pts_utm <- terra::project(pts_ll, bgi_proj)\n"
    "  bgi_vals <- terra::extract(bgi_rast, pts_utm)[,2]\n"
    "  bgi_lookup <- setNames(bgi_vals, as.character(me_plot_xy$CN))\n"
    "  bgi_match <- bgi_lookup[as.character(samp$PLT_CN_t1_str)]\n"
    "  cat(sprintf(\"[v27c] BGI extracted: n=%d matched, mean=%.3f, median=%.3f, NA=%d\\n\",\n"
    "      sum(!is.na(bgi_match)), mean(bgi_match, na.rm=TRUE),\n"
    "      median(bgi_match, na.rm=TRUE), sum(is.na(bgi_match))))\n"
    "  # Fill NAs with mean to avoid losing plots\n"
    "  bgi_mean_pos <- mean(bgi_match[bgi_match > 0], na.rm=TRUE)\n"
    "  bgi_match[is.na(bgi_match) | bgi_match <= 0] <- bgi_mean_pos\n"
    "  bgi_per_tree <- bgi_match[match(base_init$STAND, samp$PLT_CN_t1_str)]\n"
    "  bgi_per_tree[is.na(bgi_per_tree)] <- bgi_mean_pos\n"
    "} else {\n"
    "  cat(\"[v27c] terra not installed; skipping raster BGI configs\\n\")\n"
    "  bgi_per_tree <- rep(1, nrow(base_init))\n"
    "  bgi_mean_pos <- 1\n"
    "}\n"
)
s = s.replace(anchor, extract_block, 1)

# Extend the run block with two new BGI configs after fvssi_replace
import re
old_pat = re.compile(r'(stand_init <- stand_init_fvssi.*?rows\[\["c"\]\] <- summ\("fvssi_replace",   run_cfg\(TRUE, "Y", 0\)\))',
                     re.DOTALL)
m = old_pat.search(s)
if m:
    insertion = m.group(1) + (
        "\n\n"
        "# Restore CSI baseline stand_init; vary dDBH.mult by BGI per tree.\n"
        "stand_init <- stand_init_csi\n"
        "dm_orig <- base_init$dDBH.mult\n"
        "\n"
        "base_init$dDBH.mult <- bgi_per_tree\n"
        "cat(sprintf(\"[v27c] bgi_as_dmult   (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, dDBH.mult = BGI [mean %.3f])\\n\", mean(bgi_per_tree, na.rm=TRUE)))\n"
        "rows[[\"d\"]] <- summ(\"bgi_as_dmult\",   run_cfg(TRUE, \"Y\", 0))\n"
        "\n"
        "base_init$dDBH.mult <- bgi_per_tree / bgi_mean_pos\n"
        "cat(sprintf(\"[v27c] bgi_recenter   (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, dDBH.mult = BGI / %.3f)\\n\", bgi_mean_pos))\n"
        "rows[[\"e\"]] <- summ(\"bgi_recenter\",   run_cfg(TRUE, \"Y\", 0))\n"
        "\n"
        "base_init$dDBH.mult <- dm_orig"
    )
    s = s.replace(m.group(1), insertion, 1)
else:
    print("WARN: fvssi block not found, BGI configs not added")

# Update final cat header
s = s.replace(
    "FIA: do alternative site indices (SICOND, FVS_SITE_INDEX) close BA vs ClimateSI?",
    "FIA: do alternative site indices and ME BGI raster close BA vs ClimateSI?")

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "bgi_block_ok", "bgi_per_tree" in s,
      "bgi_configs_ok", "bgi_as_dmult" in s and "bgi_recenter" in s)
