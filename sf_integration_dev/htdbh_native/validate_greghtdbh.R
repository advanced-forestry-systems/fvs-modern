## validate_greghtdbh.R -- R reference for the CONUS native HT-DBH forms.
## Re-evaluates each species' exact model.formula (the six Marshall forms) over
## the same DBH grid the Fortran harness uses, then compares to htdbh_fortran.csv.
## Requires max abs diff < 1e-6 across all six forms and overall.
suppressWarnings(suppressMessages({}))
coef <- read.csv('conus_htdbh_coefficients.csv', stringsAsFactors=FALSE)
coef <- coef[coef$model_id %in% 1:6, ]
bh <- 4.5
grid <- 1:50

## the six forms, honoring the Marshall model.formula strings exactly
predHT <- function(m, B1, B2, B3, DBH) {
  if (m==1) return(bh + exp(B1 + B2*DBH^(-1.0)))
  if (m==2) return(bh + exp(B1 + B2*DBH^B3))
  if (m==3) return(bh + B1*(1.0-exp(B2*DBH))^B3)
  if (m==4) return(bh + exp(B1 + B2/(DBH+1.0)))
  if (m==5) return(bh + exp(B1 + B2/(DBH+B3)))
  if (m==6) return(bh + B1*(1.0-exp(B2*DBH^B3)))
  return(NA_real_)
}

rows <- list()
for (i in seq_len(nrow(coef))) {
  m  <- coef$model_id[i]
  B1 <- as.numeric(coef$B1[i]); B2 <- as.numeric(coef$B2[i])
  B3 <- suppressWarnings(as.numeric(coef$B3[i])); if (is.na(B3)) B3 <- 0.0
  ht <- predHT(m, B1, B2, B3, grid)
  rows[[i]] <- data.frame(SPCD=coef$SPCD[i], model_id=m, DBH=as.numeric(grid), HT_R=ht)
}
ref <- do.call(rbind, rows)
write.csv(ref, 'htdbh_reference_R.csv', row.names=FALSE)

## compare
fort <- read.csv('htdbh_fortran.csv', stringsAsFactors=FALSE)
mrg <- merge(ref, fort, by=c('SPCD','model_id','DBH'))
stopifnot(nrow(mrg) == nrow(ref))
mrg$absdiff <- abs(mrg$HT_R - mrg$HT)

cat(sprintf('rows compared: %d (species x DBH)\n', nrow(mrg)))
cat('--- max abs diff (R vs Fortran) per model form ---\n')
for (m in sort(unique(mrg$model_id))) {
  sub <- mrg[mrg$model_id==m, ]
  cat(sprintf('  model %d : n=%6d  max_abs_diff = %.3e\n', m, nrow(sub), max(sub$absdiff, na.rm=TRUE)))
}
ov <- max(mrg$absdiff, na.rm=TRUE)
cat(sprintf('--- OVERALL max abs diff = %.3e ---\n', ov))
if (is.finite(ov) && ov < 1e-6) {
  cat('PASS: max abs diff < 1e-6\n')
  quit(status=0)
} else {
  cat('FAIL: max abs diff >= 1e-6\n')
  ## show worst offenders
  print(head(mrg[order(-mrg$absdiff), c('SPCD','model_id','DBH','HT_R','HT','absdiff')], 10))
  quit(status=1)
}
