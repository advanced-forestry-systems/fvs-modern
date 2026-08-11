## Unified joint fit: one shared self-thinning response + a per-variant LEVEL scalar on the
## localized (brms) maximum SDI, estimated jointly from FIA remeasurement. The shared response
## is a logistic in relative density (identifiable against per-variant level, unlike a free
## smooth). Compares fit to (a) raw brms (all levels = 1) and (b) reports the estimated levels,
## which should track the engine-derived optimal levels.
suppressMessages({library(data.table)})
set.seed(7)
d <- as.data.table(readRDS("~/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds"))
d <- d[is.finite(SDImax_brms)&SDImax_brms>0&is.finite(SDI1)&SDI1>0&is.finite(TPH1)&is.finite(TPH2)&TPH1>0&TPH2>0&is.finite(YEARS)&YEARS>0]
d[, mort := -log(TPH2/TPH1)/YEARS]; d <- d[is.finite(mort) & mort > -0.02 & mort < 0.2]  # trim wild values
d[, region := toupper(as.character(fvs_variant))]
keep <- d[, .N, by=region][N>=1500, region]; d <- d[region %in% keep]
## subsample for optim speed, stratified by region
d <- d[, .SD[sample(.N, min(.N, 20000))], by=region]
regs <- sort(unique(d$region)); G <- length(regs)
gi <- match(d$region, regs)
SDI1 <- d$SDI1; SMX <- d$SDImax_brms; mort <- d$mort
cat("rows", nrow(d), "regions", G, "\n")
## shared logistic: mort = b0 + b1 / (1 + exp(-(RD - c)/s)), RD = SDI1/(k_g*SMX)
## params: b0,b1,c,s, then log(k_g) for each region (k>0)
negSSE <- function(p){
  b0<-p[1]; b1<-p[2]; cc<-p[3]; s<-exp(p[4]); k<-exp(p[5:(4+G)])
  RD <- SDI1/(k[gi]*SMX)
  pred <- b0 + b1/(1+exp(-(RD-cc)/s))
  sum((mort-pred)^2)
}
init <- c(0.0, 0.05, 0.6, log(0.1), rep(0, G))
fit <- optim(init, negSSE, method="BFGS", control=list(maxit=400, reltol=1e-9))
b0<-fit$par[1]; b1<-fit$par[2]; cc<-fit$par[3]; s<-exp(fit$par[4]); k<-exp(fit$par[5:(4+G)])
sse_joint <- fit$value
## baseline: same shape but all k=1 (raw brms), refit shape only
negSSE_k1 <- function(p){ b0<-p[1];b1<-p[2];cc<-p[3];s<-exp(p[4]); RD<-SDI1/SMX; pred<-b0+b1/(1+exp(-(RD-cc)/s)); sum((mort-pred)^2) }
f1 <- optim(c(0,0.05,0.6,log(0.1)), negSSE_k1, method="BFGS", control=list(maxit=400))
sst <- sum((mort-mean(mort))^2)
cat(sprintf("\nShared self-thinning curve: mort = %.4f + %.4f / (1+exp(-(RD-%.3f)/%.3f))\n", b0,b1,cc,s))
cat(sprintf("R2: raw brms (level=1) %.4f  ->  joint per-variant level %.4f\n", 1-f1$value/sst, 1-sse_joint/sst))
res <- data.table(region=regs, level_joint=round(k,2))
fwrite(res, "~/fvs-conus/output/joint_fit_levels.csv")
print(res[order(level_joint)])
cat("DONE_JOINTFIT\n")
