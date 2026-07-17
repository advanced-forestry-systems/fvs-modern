f="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v5.R"
s=open(f).read()
s=s.replace("## conus_eq_projector_v3.R -- annualized stand projection","## conus_eq_projector_v5.R (copy of v4) -- annualized stand projection")
anchor='ST_SITE_SCALE<-as.numeric(ga("st_site_scale","0.0")); SDIMAX_REF<-as.numeric(ga("sdimax_ref","935"))'
assert anchor in s,'anchor1'
s=s.replace(anchor, anchor+'\n## v5: SIZE-BIASED self-thinning. ST_BIAS>0 makes competition mortality preferentially remove smaller/low-crown trees.\nST_BIAS<-as.numeric(ga("st_bias","0.0"))')
anchor2='      s_nat<-pmin(pmax(s_nat,1e-9),1); H<- -log(s_nat)            # cycle hazard per tree\n      N1<-sum(tl$TPA)'
assert anchor2 in s,'anchor2'
inject='''      s_nat<-pmin(pmax(s_nat,1e-9),1); H<- -log(s_nat)            # cycle hazard per tree
      if(ST_BIAS>0 && length(H)>1){
        d<-tl$dbh_in; mu<-mean(d); sdv<-sd(d); if(!is.finite(sdv)||sdv<1e-6) sdv<-1
        z_size<-(d-mu)/sdv
        crown_def<-1-pmin(pmax(tl$CR,0),1)
        bias<-exp(ST_BIAS*(-z_size + crown_def))
        H<-H*bias
      }
      N1<-sum(tl$TPA)'''
s=s.replace(anchor2,inject)
anchor3='''      surv[!is.finite(surv)]<-1
      tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0'''
assert anchor3 in s,'anchor3'
diag='''      surv[!is.finite(surv)]<-1
      if(nchar(Sys.getenv("DIAG_THIN"))>0 && identical(cn,Sys.getenv("DIAG_THIN")) && cy %in% c(40,60,80,100)){
        n<-length(tl$dbh_in); ord<-order(tl$dbh_in,decreasing=TRUE)
        topn<-max(1,floor(0.10*n)); botn<-max(1,floor(0.50*n))
        topi<-ord[seq_len(topn)]; boti<-ord[seq(n-botn+1,n)]
        wmean<-function(x,w){ww<-sum(w); if(ww>0) sum(x*w)/ww else NA_real_}
        cat(sprintf("[DIAGTHIN] cn=%s cy=%d RD=%.3f | DOM(top10pct): DBH=%.1f HT=%.1f killfrac=%.4f | SUP(bot50pct): DBH=%.1f HT=%.1f killfrac=%.4f | totTPA=%.1f totKill=%.3f\\n",
          cn,cy,rd_m,
          wmean(tl$dbh_in[topi],tl$TPA[topi]),wmean(tl$HT[topi],tl$TPA[topi]),wmean((1-surv)[topi],tl$TPA[topi]),
          wmean(tl$dbh_in[boti],tl$TPA[boti]),wmean(tl$HT[boti],tl$TPA[boti]),wmean((1-surv)[boti],tl$TPA[boti]),
          sum(tl$TPA),sum(tl$TPA*(1-surv))))
      }
      tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0'''
s=s.replace(anchor3,diag)
open(f,'w').write(s)
print('patched OK')
