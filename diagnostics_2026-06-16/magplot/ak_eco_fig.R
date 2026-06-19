suppressMessages({library(data.table); library(ggplot2)})
d <- fread("/fs/scratch/PUOM0008/crsfaaron/akwork/ak_eco_validate_results.csv")
d[, dpj := t2-t0][, dob := t2o-t1o]
g <- d[, .(n=.N, incr=100*(sum(dpj)-sum(dob))/sum(dob), level=100*(sum(t2)-sum(t2o))/sum(t2o)), by=ecoL1][n>=8]
g[, ecoL1 := factor(ecoL1, levels=ecoL1[order(incr)])]
g[, analog := ifelse(grepl("MARINE",ecoL1),"AK analog (SE Alaska ecoregion)","other BC ecoregion")]
m <- melt(g, id.vars=c("ecoL1","n","analog"), measure.vars=c("incr","level"), variable.name="metric", value.name="bias")
m[, metric := ifelse(metric=="incr","BA increment bias","standing BA bias")]
p <- ggplot(m, aes(ecoL1, bias, fill=analog)) +
  geom_col(width=0.7) + geom_hline(yintercept=0, linewidth=0.3) +
  geom_text(aes(label=sprintf("%+.0f%%\n(n=%d)",bias,n)), hjust=ifelse(m$bias<0,1.05,-0.05), size=2.7) +
  facet_wrap(~metric, scales="free_x") + coord_flip() +
  scale_fill_manual(values=c("AK analog (SE Alaska ecoregion)"="#1f6f54","other BC ecoregion"="#9aa7b4")) +
  labs(title="Default FVS-Alaska bias vs observed growth, by NA Level I ecoregion (BC MAGPlot)",
       subtitle="Negative = FVS under-predicts. Clean-ingestion subset; 13-35 yr remeasurement.",
       x=NULL, y="bias (%)", fill=NULL) +
  theme_minimal(base_size=11) + theme(legend.position="top", plot.title=element_text(size=11))
png("/fs/scratch/PUOM0008/crsfaaron/akwork/fig_ak_ecoregion_bias.png", width=2400, height=1100, res=300)
print(p); dev.off()
cat("wrote figure\n")
