#!/usr/bin/env Rscript
# fig_uncertainty_layers.R -- the three uncertainty layers for Maine carbon
# density: parameter (posterior draws) vs structural (engine spread).
suppressMessages({library(data.table); library(ggplot2)})
SD <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
ci <- fread(file.path(SD, "post_ME", "posterior_ME_ci.csv"))   # calibrated param CI

# structural band: the three engines' ME reserve density (Mg C/ha), from the
# v2 aggregated series (agc_live_total reserve).
eng <- rbindlist(lapply(c("default","calibrated","gompit"), function(c){
  f <- file.path(SD, sprintf("perseus_series_%s_v2", c), "ycx_ME_state_series.csv")
  d <- fread(f)[metric=="agc_live_total" & mgmt=="reserve (no harvest)",
               .(year, value, engine=c)]
  d}), fill=TRUE)
struct <- dcast(eng, year ~ engine, value.var="value")
struct[, `:=`(lo=pmin(default,calibrated,gompit), hi=pmax(default,calibrated,gompit))]

p <- ggplot() +
  geom_ribbon(data=struct, aes(year, ymin=lo, ymax=hi, fill="structural (engine spread)"),
              alpha=0.25) +
  geom_ribbon(data=ci, aes(year, ymin=p2_5, ymax=p97_5, fill="parameter (posterior draws)"),
              alpha=0.85) +
  geom_line(data=ci, aes(year, mean), colour="#0072B2", linewidth=0.8) +
  scale_fill_manual(values=c("structural (engine spread)"="#CC9966",
                             "parameter (posterior draws)"="#0072B2")) +
  labs(x=NULL, y="live aboveground carbon (Mg C/ha)",
       title="Maine carbon: parameter vs structural uncertainty",
       subtitle="calibrated posterior-draw CI (narrow) inside the default/calibrated/gompit spread",
       fill=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", size=12))
ggsave(file.path(SD, "post_ME", "me_uncertainty_layers.png"), p,
       width=7.2, height=5, dpi=200)
cat("wrote me_uncertainty_layers.png\n")
