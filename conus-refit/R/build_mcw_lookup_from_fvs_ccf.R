#!/usr/bin/env Rscript
# Part B (MCW->CCH) FOUNDATION. Recovers per-species MAXIMUM CROWN WIDTH (MCW)
# from FVS variant CCF tables (Table 4.5.1) per David Marshall,
# "Computing Tree Maximum Crown Width Using FVS Equations" (2026-06-24).
# VERIFIED: reproduces Marshall DF example (CCF=1.03264, MCW=23.9317ft, B0=4.6389,
# B1=1.6077 @ DBH=12in, Blue Mtns DF R1=.0388 R2=.0269 R3=.00466).
# Transforms (k=0.001803026):
#  1a CCFt=R1+R2*D+R3*D^2 -> MCW=sqrt((R1+R2*D+R3*D^2)/k); linear B0=sqrt(R1/k),B1=sqrt(R3/k)
#  1b CCFt=R4*D^R5        -> A1=sqrt(R4/k), A2=R5/2; MCW=A1*D^A2
#  1c CCFt=(R1+R2+R3)*D   -> C1=sqrt((R1+R2+R3)/k); MCW=C1*sqrt(D)
# REMAINING (see report): per-variant FIA-SPCD crosswalk (DATA index->FIA), the
# CW-equation variants CA/CS/LS/NE/SN (section 4.4, cwcalc.f90), one-eq-per-species
# selection for 471 CONUS spp (genus / SW-HW fallback for ~227), then wire into
# cch_module.R and re-run ONLY the gompit arm.
suppressPackageStartupMessages({library(data.table)})
K<-0.001803026
SRC<-Sys.getenv("FVS_MODERN_SRC","/users/PUOM0008/crsfaaron/fvs-modern/src-converted")
pb<-function(L,nm){i0<-grep(paste0("^\\s*DATA\\s+",nm,"\\s*/"),L,ignore.case=TRUE);if(!length(i0))return(NULL)
 buf<-c();i<-i0[1];txt<-sub(paste0("^\\s*DATA\\s+",nm,"\\s*/"),"",L[i],ignore.case=TRUE)
 repeat{buf<-c(buf,txt);if(grepl("/",txt))break;i<-i+1;if(i>length(L))break;txt<-L[i]}
 s<-paste(buf,collapse=" ");s<-sub("/.*$","",s);s<-gsub("&"," ",s)
 v<-suppressWarnings(as.numeric(trimws(strsplit(s,",")[[1]])));v[!is.na(v)]}
ev<-function(v){f<-file.path(SRC,v,"ccfcal.f90");if(!file.exists(f))return(NULL)
 L<-readLines(f,warn=FALSE);R1<-pb(L,"RD1");R2<-pb(L,"RD2");R3<-pb(L,"RD3");R4<-pb(L,"RDA");R5<-pb(L,"RDB")
 if(is.null(R1))return(NULL);n<-length(R1)
 data.table(variant=toupper(v),sp_idx=seq_len(n),R1=R1,R2=R2[seq_len(n)],R3=R3[seq_len(n)],
  R4=if(length(R4)>=n)R4[seq_len(n)] else NA_real_,R5=if(length(R5)>=n)R5[seq_len(n)] else NA_real_,
  MCW_B0=sqrt(pmax(R1,0)/K),MCW_B1=sqrt(pmax(R3[seq_len(n)],0)/K))}
vars<-c("bm","cr","ut","ci","ec","em","ie","kt","nc","so","tt","wc","ws","ak")
out<-rbindlist(lapply(vars,ev),fill=TRUE)
of<-Sys.getenv("MCW_OUT","/users/PUOM0008/crsfaaron/fvs-conus/output/conus/mcw_ccf_coefficients_raw.csv")
dir.create(dirname(of),showWarnings=FALSE,recursive=TRUE)
fwrite(out,of)
cat(sprintf("Extracted %d species-eqs across %d CCF variants -> %s\n",nrow(out),length(unique(out[["variant"]])),of))
