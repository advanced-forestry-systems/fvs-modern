suppressPackageStartupMessages(library(data.table))
SIDIR<-'/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4'
OUT<-'/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bakuzis_singlecohort/v6_figs'
dir.create(OUT,showWarnings=FALSE)
agg_arm<-function(metdir,region){
  SI<-fread(file.path(SIDIR,paste0('standinit_',region,'.csv')),colClasses=list(character='STAND_CN'),select=c('STAND_CN','AGE','SITE_INDEX'))
  SI[,STAND_CN:=sub('\\..*$','',STAND_CN)];SI[,AGE:=as.numeric(AGE)];SI[,SITE_INDEX:=as.numeric(SITE_INDEX)];SI<-unique(SI,by='STAND_CN')
  SI[,seed:=is.finite(AGE)&AGE>=10&AGE<=40]
  M<-fread(file.path(metdir,paste0('conus_eq_',tolower(region),'_conus_b2_metrics.csv')),colClasses=list(character='STAND_CN'))
  M[,STAND_CN:=sub('\\..*$','',STAND_CN)]
  M<-merge(M,SI[,.(STAND_CN,AGE0=AGE,SI=SITE_INDEX,seed)],by='STAND_CN',all.x=TRUE)
  M<-M[seed==TRUE&is.finite(SI)&SI>0&is.finite(HT_M_DOM)]
  M[,age:=AGE0+as.integer(PROJ_YEAR)];M[,agebin:=round(age/5)*5]
  qb<-quantile(M$SI,c(0,.25,.5,.75,1),na.rm=TRUE);M[,site:=cut(SI,qb,include.lowest=TRUE,labels=1:4)]
  M[!is.na(site),.(HTDOM=mean(HT_M_DOM)),by=.(site,agebin)][order(site,agebin)]
}
png(file.path(OUT,'htdom_recession_v4_v6.png'),width=1500,height=520,res=110)
par(mfrow=c(1,3),mar=c(4,4,3,1))
for(rg in c('NE','SN','PN')){
  a4<-agg_arm('/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_v4_seed',rg)
  a6<-agg_arm('/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_v6_seed',rg)
  yl<-range(c(a4$HTDOM,a6$HTDOM),na.rm=TRUE)
  plot(NA,xlim=c(0,140),ylim=yl,xlab='cohort age (yr)',ylab='top height HT_M_DOM (m)',main=paste0(rg,': top-height (v4 solid, v6 dashed)'))
  cols<-c('1'='#762a83','2'='#1b7837','3'='#2166ac','4'='#b35806')
  for(s in 1:4){d4<-a4[site==s];d6<-a6[site==s]
    lines(d4$agebin,d4$HTDOM,col=cols[as.character(s)],lwd=2)
    lines(d6$agebin,d6$HTDOM,col=cols[as.character(s)],lwd=2,lty=2)}
  legend('bottomright',legend=paste('site',1:4),col=cols,lwd=2,bty='n',cex=0.8)
}
dev.off(); cat('htdom v4-v6 figure written\n')
