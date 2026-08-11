# Four-way comparison aggregation: variants x species x ecoregion x landowner (2026-06-18)
# Joins the per-condition default/calibrated predictions (calib_4way_*_percond.csv, keyed by CN) to the
# four keys from the remeasurement RDS, then computes default vs calibrated bias by each margin and the
# populated crosses. Output: calib_4way_margins.csv + fig_4way_landowner.png, fig_4way_ecoregion.png.
suppressWarnings(suppressMessages({library(data.table); library(ggplot2)}))
set.seed(5)
RD <- "/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16"
RDS <- "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"
elog <- function(m) cat(sprintf("[%s] %s\n", format(Sys.time()), m), file=file.path(RD,"error_log.txt"), append=TRUE)
pc_path <- Sys.glob(file.path(RD,"calib_4way_*_percond.csv"))[1]
stopifnot(!is.na(pc_path))
pc <- fread(pc_path); cat("per-condition rows:", nrow(pc), " variants:", length(unique(pc$variant)), "\n")

d <- as.data.table(readRDS(RDS))
# join key: the remeasurement RDS keys on the t1 plot (PLT_CN_cond1); the harness now emits t1.
joincol <- if ("t1" %in% names(pc)) "t1" else "cn"
cand <- intersect(c("PLT_CN_cond1","CN_cond1","PLT_CN","CN","plot_key"), names(d))
best <- NULL; bestrate <- 0
for (cc in cand) {
  r <- mean(pc[[joincol]] %in% suppressWarnings(as.numeric(d[[cc]])), na.rm=TRUE)
  cat("  match", joincol, "via", cc, ":", round(100*r,1), "%\n"); if (r > bestrate) {bestrate <- r; best <- cc}
}
if (is.null(best) || bestrate < 0.1) { elog("4way: no CN column matched > 10%"); cat("CN MATCH FAILED; columns:\n"); print(grep("CN|PLT|key", names(d), value=TRUE)); quit(status=0) }
cat("using CN column:", best, "(", round(100*bestrate,1), "% match)\n")
d[, cnj := suppressWarnings(as.numeric(get(best)))]

# per-condition keys: dominant species group, ecoregion L2, landowner group
spcol <- if ("SPGRPCD" %in% names(d)) "SPGRPCD" else "SPCD"
keys <- d[!is.na(cnj), .(
  sp = as.integer(names(sort(table(get(spcol)), decreasing=TRUE))[1]),
  ecoL2 = as.character(EPA_L2_CODE[1]),
  own = as.integer(OWNGRPCD_cond1[1])
), by = cnj]
OWNLAB <- c("10"="National Forest","20"="Other federal","30"="State/local","40"="Private")
keys[, landowner := ifelse(as.character(own) %in% names(OWNLAB), OWNLAB[as.character(own)], "Other/unknown")]

m <- merge(pc, keys, by.x=joincol, by.y="cnj", all.x=TRUE)
cat("joined rows with keys:", sum(!is.na(m$own)), "of", nrow(m), "\n")

bias <- function(pred,obs){ ok <- is.finite(pred)&is.finite(obs)&obs>0; if(sum(ok)<5) return(NA_real_); 100*sum(pred[ok]-obs[ok])/sum(obs[ok]) }
METR <- list(BA=c("dBA","cBA","oBA"), TPH=c("dTPH","cTPH","oTPH"), QMD=c("dQMD","cQMD","oQMD"), VOL=c("dVOL","cVOL","oVOL"))

# marginal comparison by each dimension
margins <- rbindlist(lapply(c("variant","sp","ecoL2","landowner"), function(dim){
  rbindlist(lapply(names(METR), function(mm){
    dd <- m[!is.na(get(dim))]
    g <- dd[, .(n=.N, def=bias(get(METR[[mm]][1]),get(METR[[mm]][3])),
                cal=bias(get(METR[[mm]][2]),get(METR[[mm]][3]))), by=dim]
    g[, `:=`(dimension=dim, metric=mm)]; setnames(g, dim, "level"); g[, level:=as.character(level)]
    g[n>=20]
  }))
}), fill=TRUE)
fwrite(margins, file.path(RD,"calib_4way_margins.csv"))
cat("wrote calib_4way_margins.csv\n")

# figure: calibrated vs default |bias| by landowner (faceted by metric)
tryCatch({
  d2 <- margins[dimension=="landowner"]
  d2l <- melt(d2, id.vars=c("level","metric","n"), measure.vars=c("def","cal"), variable.name="arm", value.name="bias")
  d2l[, arm:=ifelse(arm=="def","default","calibrated")]
  p <- ggplot(d2l, aes(level, abs(bias), fill=arm)) + geom_col(position=position_dodge(0.8), width=0.7) +
    facet_wrap(~metric, scales="free_y") + coord_flip() +
    scale_fill_manual(values=c(default="#9aa7b4", calibrated="#1f6f54")) +
    labs(title="Default vs calibrated |bias| by landowner (disturbance-clean)", x=NULL, y="|bias| (%)", fill=NULL) +
    theme_minimal(base_size=11) + theme(legend.position="top")
  png(file.path(RD,"fig_4way_landowner.png"), width=2200, height=1500, res=300); print(p); dev.off()
  png(file.path(RD,"fig_4way_landowner_thumb.png"), width=780, height=520, res=72); print(p); dev.off()
  cat("wrote fig_4way_landowner.png\n")
}, error=function(e) elog(paste("fig landowner", e$message)))

# figure: QMD calibrated bias by ecoregion x landowner heatmap (the cross)
tryCatch({
  cr <- m[!is.na(ecoL2)&!is.na(own), .(n=.N, cal_QMD=bias(cQMD,oQMD)), by=.(ecoL2, landowner)][n>=20]
  if (nrow(cr)>0) {
    p2 <- ggplot(cr, aes(landowner, ecoL2, fill=cal_QMD)) + geom_tile(colour="white") +
      geom_text(aes(label=sprintf("%+.0f",cal_QMD)), size=2.6) +
      scale_fill_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", midpoint=0) +
      labs(title="Calibrated QMD bias: ecoregion (L2) x landowner", x=NULL, y="EPA L2 ecoregion", fill="QMD bias %") +
      theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=30,hjust=1))
    png(file.path(RD,"fig_4way_ecoregion.png"), width=2200, height=1700, res=300); print(p2); dev.off()
    png(file.path(RD,"fig_4way_ecoregion_thumb.png"), width=780, height=600, res=72); print(p2); dev.off()
    cat("wrote fig_4way_ecoregion.png\n")
  }
}, error=function(e) elog(paste("fig eco", e$message)))
gc(); cat("DONE_4WAY_AGG\n")
