## debug_ingrowth.R - trace exactly where ingrowth dies (CutPoint=0 expected value)
suppressWarnings(suppressMessages({ library(plyr); library(dplyr); library(purrr); library(tibble) }))
source(Sys.getenv("ACD126","/users/PUOM0008/crsfaaron/AcadianGY_12.3.6.r"))
cat("model:", AcadianVersionTag, "\n")

set.seed(7); n <- 60
sp <- sample(c("RS","BF","RM","YB","WP","SM"), n, replace=TRUE, prob=c(.28,.22,.20,.12,.08,.10))
tr <- data.frame(PLOT=1L, TREE=seq_len(n), SP=sp, DBH=pmax(3, rgamma(n,4,0.22)),
  CR=pmin(.9,pmax(.2,rnorm(n,.55,.12))), EXPF=runif(n,3,12), stringsAsFactors=FALSE)
tr$HT <- 1.37 + 0.85*tr$DBH^0.78; tr$HCB <- tr$HT*(1-tr$CR)
tr$max.dbh <- 200; tr$max.height <- 60; tr$dDBH.mult <- 1; tr$dHt.mult <- 1; tr$mort.mult <- 1
ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=3.0, CutPoint=0)
st <- list(CSI=14, ELEV=200)

tr$YEAR <- 1
o1 <- AcadianGYOneStand(tr, stand=st, ops=ops)
cat(sprintf("\ncycle1: in=%d out=%d  (added %d)\n", nrow(tr), nrow(o1), nrow(o1)-nrow(tr)))
cat("out1 columns:", paste(names(o1), collapse=","), "\n")
newsm <- o1[o1$DBH <= 4 & !is.na(o1$DBH), ]
cat(sprintf("out1 small trees DBH<=4: %d ; min DBH=%.2f\n", nrow(newsm), min(o1$DBH, na.rm=TRUE)))
for (c0 in c("dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height"))
  if (c0 %in% names(o1)) cat(sprintf("  out1 %-10s NA count = %d\n", c0, sum(is.na(o1[[c0]]))))

o1$YEAR <- 2
o2 <- AcadianGYOneStand(o1, stand=st, ops=ops)
cat(sprintf("\ncycle2: in=%d out=%d  NA-DBH in out2=%d\n", nrow(o1), nrow(o2), sum(is.na(o2$DBH))))
ba <- function(d) sum(d$EXPF*0.00007854*d$DBH^2, na.rm=TRUE)
cat(sprintf("BA: start=%.2f  c1=%.2f  c2=%.2f\n", ba(tr), ba(o1), ba(o2)))
