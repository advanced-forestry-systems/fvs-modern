import pandas as pd, numpy as np
f="/users/PUOM0008/crsfaaron/fvs-conus/output/comparisons/intermediate/validation_data.csv.gz"
d=pd.read_csv(f, compression="gzip", low_memory=False)
def analyze(df,obs,pred,lo,hi,label):
    m=df[[obs,pred,lo,hi]].dropna()
    m=m[m[obs]>0]
    if len(m)<30: return None
    o=m[obs].values; p=m[pred].values
    cur_cov=((o>=m[lo].values)&(o<=m[hi].values)).mean()*100
    cur_hw=((m[hi].values-m[lo].values)/2)
    resid=o-p
    sd=resid.std()
    # corrected predictive interval: pred +/- 1.96*residual SD
    nlo=p-1.96*sd; nhi=p+1.96*sd
    new_cov=((o>=nlo)&(o<=nhi)).mean()*100
    return dict(n=len(m), cur_cov=cur_cov, cur_hw_pct=100*np.mean(cur_hw)/o.mean(),
                resid_sd_pct=100*sd/o.mean(), needed_hw_pct=100*1.96*sd/o.mean(), new_cov=new_cov)
print("BA: current vs residual-augmented predictive interval coverage")
print("%-9s %6s %9s %10s %11s %9s" % ("variant","n","cur_PICP","cur_HW%","resid_SD%","new_PICP"))
for v in ["OVERALL"]+sorted(d.VARIANT.dropna().unique().tolist()):
    sub=d if v=="OVERALL" else d[d.VARIANT==v]
    r=analyze(sub,"BA_t2","BA_pred_calib","BA_pred_calib_lo","BA_pred_calib_hi","BA")
    if r: print("%-9s %6d %8.1f%% %9.1f %10.1f %8.1f%%" % (v,r["n"],r["cur_cov"],r["cur_hw_pct"],r["resid_sd_pct"],r["new_cov"]))
