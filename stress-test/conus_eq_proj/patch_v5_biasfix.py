f="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v5.R"
s=open(f).read()
old='''      if(ST_BIAS>0 && length(H)>1){
        d<-tl$dbh_in; mu<-mean(d); sdv<-sd(d); if(!is.finite(sdv)||sdv<1e-6) sdv<-1
        z_size<-(d-mu)/sdv
        crown_def<-1-pmin(pmax(tl$CR,0),1)
        bias<-exp(ST_BIAS*(-z_size + crown_def))
        H<-H*bias
      }'''
new='''      if(ST_BIAS>0 && length(H)>1){
        d<-tl$dbh_in; mu<-mean(d); sdv<-sd(d); if(!is.finite(sdv)||sdv<1e-6) sdv<-1
        z_size<-pmin(pmax((d-mu)/sdv,-3),3)
        crown_def<-1-pmin(pmax(tl$CR,0),1)
        bias<-exp(pmin(pmax(ST_BIAS*(-z_size + crown_def),-3),3))   # bounded multiplier
        bias[!is.finite(bias)]<-1
        H<-H*bias; H<-pmin(pmax(H,0),50)
      }'''
assert old in s, 'biasfix anchor missing'
s=s.replace(old,new)
open(f,'w').write(s)
print('biasfix applied')
