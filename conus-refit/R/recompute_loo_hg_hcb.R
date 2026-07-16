suppressMessages({library(loo)})
comp <- function(label, dir, bname, tname, base_lab, test_lab){
  bf <- file.path(dir, paste0(bname,"_fit.rds")); tf <- file.path(dir, paste0(tname,"_fit.rds"))
  if(!file.exists(bf)||!file.exists(tf)){cat(label,": missing fit(s)\n");return(invisible())}
  cat("==",label,"== loading base\n")
  fb <- readRDS(bf); lb <- loo(fb$draws("log_lik", format="draws_matrix")); rm(fb); gc()
  saveRDS(lb, file.path(dir, paste0(bname,"_loo.rds")))
  cat("== loading test\n")
  ft <- readRDS(tf); lt <- loo(ft$draws("log_lik", format="draws_matrix")); rm(ft); gc()
  saveRDS(lt, file.path(dir, paste0(tname,"_loo.rds")))
  cmp <- loo_compare(setNames(list(lb,lt), c(base_lab,test_lab)))
  cat("\n---- ",label," loo_compare ----\n"); print(cmp)
  se <- round(abs(cmp[2,1]/cmp[2,2]),2)
  cat("\nSEs:", se, "  winner:", rownames(cmp)[1], "\n\n")
}
comp("HG", "output/conus/hg/v6testF", "hg_b","hg_t","bgi","plus_v6")
cat("DONE_HG_FAST\n")
comp("HCB","output/conus/hcb/v6testF","hcb_b","hcb_t","cspiv4","plus_v6")
cat("DONE_HCB_FAST\n")
