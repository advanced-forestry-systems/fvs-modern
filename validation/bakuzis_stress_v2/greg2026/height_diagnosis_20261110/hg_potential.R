suppressPackageStartupMessages(library(data.table))
hg <- fread("/users/PUOM0008/crsfaaron/fvs-modern-wt-greg2026/config/greg_hg_coefficients.csv")
cl <- fread("/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/climate/maine_greg_pc_climate_summary_1991-2020.csv")[stat=="median"]
sp <- data.table(code=c("SM","YB","AB","WA","RO","WP","RM","WO","BF","RS","BS"), SPCD=c(318,371,531,541,833,129,316,802,12,97,95),
                 grp=c(rep("Northern-Hardwood",4),"Oak-Pine","Oak-Pine","Oak-Pine","Oak-Pine","Spruce-Fir","Spruce-Fir","Spruce-Fir"))
hgf <- function(h, ht, cr, ccfl, cch, pc1, pc2) with(h, B0*B1*B2*cr^B3*exp(-B1*ht - B4*ccfl - B5*sqrt(cch) + B6*pc1 + B7*pc2)*(1-exp(-B1*ht))^(B2-1))
res <- rbindlist(lapply(seq_len(nrow(sp)), function(i) {
  h <- hg[SPCD == sp$SPCD[i]]; c0 <- cl[group == sp$grp[i]]
  rbindlist(lapply(list(c(0.6,0,0.2), c(0.45,40,0.4)), function(st) {
    ht <- 25; for (y in 1:100) ht <- ht + hgf(h, ht, st[1], st[2], st[3], c0$PC1, c0$PC2)
    data.table(sp$code[i], sp$grp[i], B0 = h$B0, state = sprintf("cr%.2f_ccfl%d", st[1], st[2]), ht100 = round(ht, 1))
  }))
}))
fwrite(res, "/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/ooda_20261110/hg_potential_100yr.csv"); print(dcast(res, V1+V2+B0~state, value.var="ht100"))
