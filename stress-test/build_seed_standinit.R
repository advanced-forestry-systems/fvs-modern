#!/usr/bin/env Rscript
# Build a seed-cohort standinit subset (AGE in [10,40]) for NE/SN/PN, stratified
# by within-region SITE_INDEX quartile, capped per region. Keeps full standinit
# columns so the engine driver runs unchanged. Writes a manifest of (idx,variant,batch,batchsize).
suppressPackageStartupMessages(library(data.table))
SCR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
SIDIR<-file.path(SCR,"standinit_by_variant")
OUTDIR<-file.path(SCR,"standinit_seed_v4"); dir.create(OUTDIR,showWarnings=FALSE)
CAP<-1800L; BATCH<-600L; set.seed(7)
regions<-c("NE","SN","PN")
man<-list(); idx<-0L
for(v in regions){
  d<-fread(file.path(SIDIR,paste0("standinit_",v,".csv")), colClasses=list(character="STAND_CN"))
  d[, AGEn := suppressWarnings(as.numeric(AGE))]
  d[, SIn  := suppressWarnings(as.numeric(SITE_INDEX))]
  seed<-d[is.finite(AGEn)&AGEn>=10&AGEn<=40 & is.finite(SIn)&SIn>0]
  qb<-quantile(seed$SIn, c(0,.25,.5,.75,1), na.rm=TRUE)
  seed[, sc := cut(SIn, breaks=qb, include.lowest=TRUE, labels=1:4)]
  per<-ceiling(CAP/4)
  samp<-seed[, .SD[sample(.N, min(.N, per))], by=sc]
  samp[, c("AGEn","SIn","sc"):=NULL]
  fwrite(samp, file.path(OUTDIR, paste0("standinit_",v,".csv")))
  nb<-ceiling(nrow(samp)/BATCH)
  for(b in 0:(nb-1)){ idx<-idx+1L; man[[length(man)+1]]<-data.table(idx=idx, variant=v, batch=b, batchsize=BATCH) }
  cat(sprintf("%s: seed-sample=%d  batches=%d\n", v, nrow(samp), nb))
}
M<-rbindlist(man)
fwrite(M, file.path(SCR,"seed_v4_manifest.tsv"), sep="\t", col.names=FALSE)
cat("total tasks:", nrow(M), "\n")
