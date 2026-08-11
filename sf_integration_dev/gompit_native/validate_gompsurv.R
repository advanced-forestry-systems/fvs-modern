#!/usr/bin/env Rscript
# Compare the engine GOMPSURV output (surv_fortran.csv) to the identical gompit
# survival computed in R from the same coefficient CSV. Validates that the CSV
# parses correctly into b0..b4 and that the native formula matches the reference.
co <- read.csv("greg_mortality_coefficients.csv")
fo <- read.csv("surv_fortran.csv")
g <- function(b0,b1,b2,b3,b4,cr,cch,fint){
  crc <- pmin(pmax(cr,1e-4),1); cchc <- pmax(cch,0)
  cterm <- ifelse(cchc>0, cchc^b4, 0)
  eta <- b0 + b1*(crc+0.01)^b2 + b3*cterm
  eta <- pmin(pmax(eta,-30),30)
  exp(-exp(eta)*fint)
}
m <- merge(fo, co, by="SPCD")
m$surv_r <- with(m, g(b0,b1,b2,b3,b4,cr,cch,fint))
m$absdiff <- abs(m$surv - m$surv_r)
cat(sprintf("rows compared: %d  | species: %d\n", nrow(m), length(unique(m$SPCD))))
cat(sprintf("max abs diff (Fortran single vs R double): %.3e   mean: %.3e\n",
            max(m$absdiff), mean(m$absdiff)))
cat(sprintf("survival range: [%.4f, %.4f]\n", min(m$surv_r), max(m$surv_r)))
ok <- max(m$absdiff) < 1e-4
cat(if (ok) "PASS: native GOMPSURV reproduces the R gompit to single-precision\n"
    else "FAIL: discrepancy exceeds 1e-4 -- investigate\n")
print(head(m[order(-m$absdiff), c("SPCD","cr","cch","surv","surv_r","absdiff")], 4))
