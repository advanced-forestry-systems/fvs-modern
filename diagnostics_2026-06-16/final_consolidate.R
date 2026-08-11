# Finalization consolidation for fvs-modern x fvs-conus (2026-06-18)
# Reads the committed result CSVs, builds the master summary tables, and renders three headless figures.
# Inputs (diagnostics_2026-06-16/): fourarm_abcd_20260618.csv, fourarm_engine_20260618.csv,
#   fourarm_projector_NE_20260618.csv, held_out_density_dependent_20260618.csv,
#   brms_match_rate_20260618.csv, mcw_cross_variant_spread.csv, mcw_conus_consensus.csv
# Outputs: final_master_byarm.csv, final_oos_qmd.csv, fig_fourarm_byarm.png, fig_qmd_oos.png,
#   fig_mcw_spread.png  (PNG, 300 dpi, headless). All IO wrapped; failures append to error_log.txt.
suppressWarnings(suppressMessages({library(data.table); library(ggplot2)}))
set.seed(5)
DIR <- Sys.getenv("FINDIR", "/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16")
OUT <- Sys.getenv("FINOUT", DIR)
elog <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time()), msg), file=file.path(OUT,"error_log.txt"), append=TRUE)
rd <- function(f) tryCatch(fread(file.path(DIR,f)), error=function(e){elog(paste("read fail",f,e$message)); NULL})
theme_set(theme_minimal(base_size=12))
PAL <- c(A="#9aa7b4", B="#1f6f54", C="#c47f3d", D="#2e5f8a")

# ---- load ----
abcd <- rd("fourarm_abcd_20260618.csv")
eng  <- rd("fourarm_engine_20260618.csv")
proj <- rd("fourarm_projector_NE_20260618.csv")
ho   <- rd("held_out_density_dependent_20260618.csv")
brms <- rd("brms_match_rate_20260618.csv")
mcw  <- rd("mcw_cross_variant_spread.csv")

ci_lo <- function(s) as.numeric(sub("\\[([-+0-9.]+),.*","\\1", s))
ci_hi <- function(s) as.numeric(sub(".*,([-+0-9.]+)\\]","\\1", s))

# ---- master median |bias| by arm x metric (primary = ABCD, 21 variants) ----
master <- rbindlist(lapply(c("BA","TPH","QMD","VOL"), function(m){
  cols <- paste0(c("A_","B_","C_","D_"), m)
  cols <- cols[cols %in% names(abcd)]
  if(!length(cols)) return(NULL)
  data.table(metric=m,
    arm=sub("_.*","",cols),
    median_abs_bias=sapply(cols, function(cc) median(abs(abcd[[cc]]), na.rm=TRUE)),
    n_variants=sapply(cols, function(cc) sum(is.finite(abcd[[cc]]))))
}), fill=TRUE)
fwrite(master, file.path(OUT,"final_master_byarm.csv"))
cat("MASTER median |bias| by arm x metric (ABCD, n variants per cell):\n"); print(master)

# ---- per-variant OOS QMD: default (A) vs calibrated (B) with bootstrap CIs (engine, 8 variants) ----
if(!is.null(eng)){
  oos <- eng[, .(variant, A_QMD, B_QMD,
                 B_lo=ci_lo(B_QMD_ci), B_hi=ci_hi(B_QMD_ci),
                 A_lo=ci_lo(A_QMD_ci), A_hi=ci_hi(A_QMD_ci))]
  fwrite(oos, file.path(OUT,"final_oos_qmd.csv"))
} else oos <- NULL

# ---- figure 1: median |bias| by metric x arm ----
tryCatch({
  d <- copy(master); d[, arm:=factor(arm, levels=c("A","B","C","D"))]
  d[, metric:=factor(metric, levels=c("BA","TPH","QMD","VOL"))]
  p1 <- ggplot(d, aes(metric, median_abs_bias, fill=arm)) +
    geom_col(position=position_dodge(0.8), width=0.75, colour="white") +
    geom_text(aes(label=sprintf("%.1f",median_abs_bias)),
              position=position_dodge(0.8), vjust=-0.3, size=3) +
    scale_fill_manual(values=PAL,
      labels=c(A="A default", B="B density layer", C="C growth emul.", D="D combined")) +
    labs(title="Four-arm median |bias|, FVS engine, out-of-sample (21 variants)",
         subtitle="Density layer (B) is the workhorse; growth emulation (C) adds no median gain; D tracks B",
         x="metric", y="median |bias|  (%)", fill=NULL) +
    theme(legend.position="top")
  png(file.path(OUT,"fig_fourarm_byarm.png"), width=2400, height=1500, res=300); print(p1); dev.off()
  png(file.path(OUT,"fig_fourarm_byarm_thumb.png"), width=800, height=500, res=72); print(p1); dev.off()
  cat("wrote fig_fourarm_byarm.png\n")
}, error=function(e) elog(paste("fig1 fail", e$message)))

# ---- figure 2: per-variant QMD bias, default vs calibrated, with CIs (OOS) ----
tryCatch({
  if(!is.null(oos)){
    m <- melt(oos, id.vars="variant", measure.vars=c("A_QMD","B_QMD"),
              variable.name="arm", value.name="bias")
    m[, arm:=ifelse(arm=="A_QMD","A default","B calibrated")]
    m[oos, on="variant", `:=`(lo=ifelse(arm=="B calibrated", i.B_lo, i.A_lo),
                              hi=ifelse(arm=="B calibrated", i.B_hi, i.A_hi))]
    m[, variant:=factor(variant, levels=oos[order(B_QMD)]$variant)]
    p2 <- ggplot(m, aes(variant, bias, colour=arm)) +
      geom_hline(yintercept=0, linewidth=0.3, colour="grey60") +
      geom_errorbar(aes(ymin=lo, ymax=hi), width=0.25, position=position_dodge(0.5)) +
      geom_point(size=2.2, position=position_dodge(0.5)) +
      scale_colour_manual(values=c("A default"="#9aa7b4","B calibrated"="#1f6f54")) +
      labs(title="QMD bias by variant, out-of-sample (held-out spatial fold)",
           subtitle="Calibration (brms maxSDI + density recruitment + BAIMULT) collapses QMD bias toward zero; bars are 95% bootstrap CIs",
           x=NULL, y="QMD bias  (%)", colour=NULL) +
      theme(legend.position="top")
    png(file.path(OUT,"fig_qmd_oos.png"), width=2200, height=1400, res=300); print(p2); dev.off()
    png(file.path(OUT,"fig_qmd_oos_thumb.png"), width=760, height=480, res=72); print(p2); dev.off()
    cat("wrote fig_qmd_oos.png\n")
  }
}, error=function(e) elog(paste("fig2 fail", e$message)))

# ---- figure 3: cross-variant MCW spread by species (top spread) ----
tryCatch({
  if(!is.null(mcw)){
    d <- mcw[n_variants>1][order(-MCW20_range)][1:min(15,.N)]
    d[, sp_abbrev:=factor(sp_abbrev, levels=rev(sp_abbrev))]
    p3 <- ggplot(d, aes(y=sp_abbrev)) +
      geom_segment(aes(x=MCW20_min, xend=MCW20_max, yend=sp_abbrev), colour="#c47f3d", linewidth=1.4) +
      geom_point(aes(x=MCW20_mean), colour="#1f3a5f", size=2.3) +
      labs(title="Cross-variant maximum crown width spread at 20 in DBH",
           subtitle="Same species, different open-grown crown width across western CCF variants; dot is the mean (consensus target)",
           x="MCW at 20 in DBH  (ft)", y=NULL)
    png(file.path(OUT,"fig_mcw_spread.png"), width=2200, height=1400, res=300); print(p3); dev.off()
    png(file.path(OUT,"fig_mcw_spread_thumb.png"), width=760, height=480, res=72); print(p3); dev.off()
    cat("wrote fig_mcw_spread.png\n")
  }
}, error=function(e) elog(paste("fig3 fail", e$message)))

gc()
cat("DONE_CONSOLIDATE\n")
