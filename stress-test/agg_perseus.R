suppressMessages({library(readr); library(dplyr)})
D <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_perseus_wo1"
fs <- list.files(D, pattern="^perseus_100yr_agb_batch.*\\.csv$", full.names=TRUE)
cat("files:", length(fs), "\n")
dat <- bind_rows(lapply(fs, function(f) suppressMessages(read_csv(f, show_col_types=FALSE))))
cat("rows:", nrow(dat), "\n")
write_csv(dat, file.path(D, "perseus_100yr_agb_all.csv"))
summ <- dat %>% group_by(YEAR, VARIANT, CONFIG) %>%
  summarise(n_plots=n_distinct(PLOT), mean_agb_tons_ac=mean(AGB_TONS_AC, na.rm=TRUE), .groups="drop") %>%
  arrange(VARIANT, CONFIG, YEAR)
write_csv(summ, file.path(D, "perseus_100yr_agb_summary.csv"))
cat("=== summary rows:", nrow(summ), "===\n")
print(as.data.frame(summ), row.names=FALSE)
