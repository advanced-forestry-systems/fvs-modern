#!/usr/bin/env Rscript
# Gompit-in-FVS across eight variants (four mortality-routine families):
#   eastern (vls/morts.f90): NE, CS, LS    eastern (own): SN
#   western (own Dixon/VARMRT): CR, WS, EC, CA
# Reads the per-variant native/gompit validation CSVs and produces two panels:
#   A  small-multiple AGB trajectories (native vs gompit), free y per variant
#   B  percent change in AGB at projection year 100, by variant, coloured by
#      region -- the headline that the in-engine substitution is bounded and
#      behaves variant-sensibly everywhere it was wired.
#
# Usage: Rscript 42_gompit_fvs_allvariants_figure.R --indir <dir> --out <png>

suppressWarnings(suppressMessages({library(tidyverse); library(patchwork)}))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f,d=NULL){i<-grep(paste0("^--",f,"="),args,value=TRUE);if(!length(i))return(d);sub(paste0("^--",f,"="),"",i[1])}
INDIR <- ga("indir","."); OUT <- ga("out","gompit_fvs_allvariants.png")

vorder <- c("ne","cs","ls","sn","cr","ws","ec","ca")
region <- c(ne="East",cs="East",ls="East",sn="East",cr="West",ws="West",ec="West",ca="West")
vlab <- c(ne="NE",cs="CS",ls="LS",sn="SN",cr="CR",ws="WS",ec="EC",ca="CA")

files <- list.files(INDIR, pattern="^val_(ne|cs|ls|sn|cr|ws|ec|ca)_(native|gompit)\\.csv$", full.names=TRUE)
dat <- map_dfr(files, read_csv, show_col_types=FALSE) %>%
  mutate(v=tolower(VARIANT),
         variant=factor(v, levels=vorder, labels=vlab[vorder]),
         mode=factor(MODE, levels=c("native","gompit"), labels=c("Native FVS","Gompit-in-FVS")))

summ <- dat %>% group_by(variant, mode, PROJ_YEAR) %>%
  summarise(AGB=mean(AGB), .groups="drop")

col_nat<-"#999999"; col_gom<-"#0072B2"
pA <- ggplot(summ, aes(PROJ_YEAR, AGB, colour=mode)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~variant, nrow=2, scales="free_y") +
  scale_colour_manual(values=c("Native FVS"=col_nat,"Gompit-in-FVS"=col_gom), name=NULL) +
  labs(x="Projection year", y="Mean AGB (tons/ac)",
       title="A  In-engine gompit vs native, eight variants",
       subtitle="100-yr no-harvest; bounded and realistic everywhere, no runaway") +
  theme_minimal(base_size=11) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"), strip.text=element_text(face="bold"))

chg <- dat %>% filter(PROJ_YEAR==100) %>%
  group_by(v, variant, mode) %>% summarise(AGB=mean(AGB), .groups="drop") %>%
  pivot_wider(names_from=mode, values_from=AGB) %>%
  mutate(pct=(`Gompit-in-FVS`/`Native FVS`-1)*100,
         Region=region[v]) %>%
  arrange(pct)
chg$variant <- factor(chg$variant, levels=chg$variant)

pB <- ggplot(chg, aes(variant, pct, fill=Region)) +
  geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%+.0f%%", pct)), vjust=ifelse(chg$pct<0,1.2,-0.4), size=3.2) +
  scale_fill_manual(values=c(East="#D55E00", West="#009E73")) +
  scale_y_continuous(expand=expansion(mult=c(0.12,0.08))) +
  labs(x=NULL, y="AGB change at yr 100 (gompit vs native)",
       title="B  Variant-specific effect",
       subtitle="EC (PNW conifer) smallest: ORGANON crown proxy best-fit there") +
  theme_minimal(base_size=11) +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"), legend.position="top")

fig <- pA / pB + plot_layout(heights=c(1.4,1))
ggsave(OUT, fig, width=10, height=8.4, dpi=300, bg="white")
cat("wrote", OUT, "\n"); print(chg %>% select(variant, Region, pct) %>% mutate(pct=round(pct)))
