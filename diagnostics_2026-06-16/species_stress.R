# fvs-conus species-specific vs species-free stress test (2026-06-18)
# For each component, quantify the species-specific contribution: sigma_sp (species random-intercept SD on
# the linear-predictor scale), its response-scale multiplier exp(sigma_sp), whether the species term is
# well identified (mean/sd of sigma_sp), and the fraction of species whose intercept is credibly nonzero.
# Inputs: <comp>_cspi_traits1_*_species_intercepts.csv and *_fixed_summary.csv in output/conus/.
# Output: species_stress_master.csv + fig_species_stress.png (headless, 300 dpi). Errors -> error_log.txt.
suppressWarnings(suppressMessages({library(data.table); library(ggplot2)}))
set.seed(5)
D <- Sys.getenv("CONUSDIR", "/users/PUOM0008/crsfaaron/fvs-conus/output/conus")
OUT <- Sys.getenv("FINOUT", "/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16")
elog <- function(m) cat(sprintf("[%s] %s\n", format(Sys.time()), m), file=file.path(OUT,"error_log.txt"), append=TRUE)
rd <- function(p) tryCatch(fread(p), error=function(e){elog(paste("read",p,e$message)); NULL})

# component: prefix, link scale, label
COMP <- list(
  dg   = list(pre="dg_kuehne_cspi_traits1",        link="log",   label="Diameter growth"),
  hg   = list(pre="hg_organon_fixedK_cspi_traits1", link="log",   label="Height growth"),
  hcb  = list(pre="hcb_organon_cspi_traits1",       link="log",   label="Height to crown base"),
  cr   = list(pre="cr_recession_cspi_traits1",      link="log",   label="Crown recession"),
  htdbh= list(pre="htdbh_wykoff_lognormal_cspi_traits1", link="log", label="Height-diameter"),
  mort = list(pre="mort_logit_simple_cspi_traits1", link="logit", label="Mortality")
)
get1 <- function(dt, v) { if(is.null(dt)) return(NA_real_); r<-dt[variable==v]; if(nrow(r)) as.numeric(r$mean[1]) else NA_real_ }
get1sd<- function(dt, v) { if(is.null(dt)) return(NA_real_); r<-dt[variable==v]; if(nrow(r)) as.numeric(r$sd[1]) else NA_real_ }

rows <- rbindlist(lapply(names(COMP), function(k){
  cc <- COMP[[k]]
  fx <- rd(file.path(D, paste0(cc$pre, "_fixed_summary.csv")))
  si <- rd(file.path(D, paste0(cc$pre, "_species_intercepts.csv")))
  sigma_sp <- get1(fx,"sigma_sp"); sigma_sp_sd <- get1sd(fx,"sigma_sp")
  resid <- get1(fx,"sigma"); if(is.na(resid)) resid <- get1(fx,"phi")
  # realized per-species spread from the intercept draws (z_sp are standardized; scale by sigma_sp)
  nsp <- if(!is.null(si)) nrow(si) else NA_integer_
  frac_nz <- if(!is.null(si) && all(c("mean","sd") %in% names(si))) mean(abs(si$mean)/si$sd > 2, na.rm=TRUE) else NA_real_
  realized_sd <- if(!is.null(si)) sd(si$mean, na.rm=TRUE) * ifelse(is.na(sigma_sp),1,sigma_sp) else NA_real_
  # response-scale interpretation: exp(sigma_sp) is the typical multiplicative per-species deviation
  mult <- if(!is.na(sigma_sp)) exp(sigma_sp) else NA_real_
  identified <- if(!is.na(sigma_sp) && !is.na(sigma_sp_sd) && sigma_sp_sd>0) sigma_sp/sigma_sp_sd else NA_real_
  data.table(component=k, label=cc$label, link=cc$link, n_species=nsp,
             sigma_sp=round(sigma_sp,3), sigma_sp_sd=round(sigma_sp_sd,3),
             species_mult=round(mult,3), pct_dev=round(100*(mult-1),1),
             identified_ratio=round(identified,1), frac_species_nonzero=round(frac_nz,2),
             residual_scale=round(resid,3))
}), fill=TRUE)
setorder(rows, -sigma_sp)
fwrite(rows, file.path(OUT,"species_stress_master.csv"))
cat("SPECIES-SPECIFIC contribution by component (sigma_sp on linear predictor; species_mult = exp):\n"); print(rows)

# figure: response-scale per-species deviation by component, flagged by identifiability
tryCatch({
  d <- copy(rows)[!is.na(pct_dev)]
  d[, label:=factor(label, levels=label[order(pct_dev)])]
  d[, well_id:=ifelse(identified_ratio>=2,"identified (mean/sd>=2)","weakly identified")]
  p <- ggplot(d, aes(label, pct_dev, fill=well_id)) +
    geom_col(width=0.7, colour="white") +
    geom_text(aes(label=sprintf("%.0f%%", pct_dev)), hjust=-0.15, size=3.4) +
    coord_flip() +
    scale_fill_manual(values=c("identified (mean/sd>=2)"="#1f6f54","weakly identified"="#c47f3d")) +
    labs(title="fvs-conus: typical per-species deviation captured by the species-specific term",
         subtitle="exp(sigma_sp) - 1 on the response scale (odds for mortality). Larger = more lost by going species-free.",
         x=NULL, y="typical per-species deviation (%)", fill=NULL) +
    expand_limits(y=max(d$pct_dev)*1.18) + theme_minimal(base_size=12) + theme(legend.position="top")
  png(file.path(OUT,"fig_species_stress.png"), width=2200, height=1300, res=300); print(p); dev.off()
  png(file.path(OUT,"fig_species_stress_thumb.png"), width=780, height=460, res=72); print(p); dev.off()
  cat("wrote fig_species_stress.png\n")
}, error=function(e) elog(paste("fig fail", e$message)))
gc(); cat("DONE_SPECIES_STRESS\n")
