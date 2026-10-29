#!/usr/bin/env Rscript
## prep_scorecard_long.R (2026-10-29)
## Attach site and density classes to the catalog engine trajectories so the v2 scorer can evaluate the site-ordering,
## Eichhorn, Sukachev and yield-density relations per variant. Classes are within-variant terciles of the condition's
## FIA site index (SICOND) and live basal area (BALIVE) among seed-passing stands, read from the seed gate file.
## Projector arms are appended unchanged (site "unspecified", density "obs"): those relations stay NOT EVALUATED there.
suppressMessages(library(data.table))
args <- commandArgs(TRUE)
ga <- function(k, d = NA) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v else d }
L <- fread(ga("long")); S <- fread(ga("seedgate"), colClasses = list(character = "plt_cn")); OUT <- ga("out", "catalog_long_scorecard.csv")
L[, plt_cn := as.character(plt_cn)]
L <- L[seed_ok %in% c(TRUE, "TRUE", "True")]
S <- S[status == "OK" & config == "default" & arm == "native", .(group = toupper(variant), plt_cn, sicond = as.numeric(sicond), balive = as.numeric(balive), stdage = as.numeric(stdage))]
S <- unique(S, by = c("group", "plt_cn"))
tercile <- function(x) { q <- quantile(x, c(1/3, 2/3), na.rm = TRUE); ifelse(is.na(x), NA, ifelse(x <= q[1], "Low", ifelse(x <= q[2], "Medium", "High"))) }
S[, `:=`(site = tercile(sicond), density = tercile(balive)), by = group]
cat("stands with classes:", nrow(S), " site NA:", sum(is.na(S$site)), " density NA:", sum(is.na(S$density)), "\n")
L <- merge(L[, !c("site", "density"), with = FALSE], S[, .(group, plt_cn, site, density)], by = c("group", "plt_cn"), all.x = TRUE)
L[is.na(site), site := "unspecified"]; L[is.na(density), density := "obs"]
proj <- rbindlist(lapply(strsplit(ga("proj"), ",")[[1]], fread), fill = TRUE)
need <- c("arm", "group", "site", "density", "traj", "age", "TPH", "BAPH", "QMD", "HT", "VOL")
out <- rbind(L[, ..need], proj[, ..need], fill = TRUE)
fwrite(out, OUT)
cat("wrote", OUT, nrow(out), "rows;", uniqueN(out$group), "groups;", uniqueN(out$arm), "arms\n")
print(dcast(unique(L[, .(group, traj, site)])[, .N, by = .(group, site)], group ~ site, value.var = "N"))
