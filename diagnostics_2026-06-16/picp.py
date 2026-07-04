import pandas as pd, numpy as np
f="/users/PUOM0008/crsfaaron/fvs-conus/output/comparisons/intermediate/validation_data.csv.gz"
d=pd.read_csv(f, compression="gzip", low_memory=False)
print("total conditions:", len(d))
def picp(df,obs,lo,hi):
    m=df[[obs,lo,hi]].dropna()
    if len(m)==0: return (float("nan"),0)
    cov=((m[obs]>=m[lo])&(m[obs]<=m[hi])).mean()
    return (100*cov, len(m))
print("\n=== interval coverage (PICP), calibrated ===")
print("%-9s %9s %10s %9s %8s" % ("variant","BA_PICP%","VOL_PICP%","HT_PICP%","n"))
for v in ["OVERALL"]+sorted(d.VARIANT.dropna().unique().tolist()):
    sub=d if v=="OVERALL" else d[d.VARIANT==v]
    ba=picp(sub,"BA_t2","BA_pred_calib_lo","BA_pred_calib_hi")
    vol=picp(sub,"VOL_CFGRS_t2","VOL_CFGRS_pred_calib_lo","VOL_CFGRS_pred_calib_hi")
    ht=picp(sub,"HT_top_t2","HT_top_calib_lo","HT_top_calib_hi")
    print("%-9s %9.1f %10.1f %9.1f %8d" % (v,ba[0],vol[0],ht[0],ba[1]))
m=d[["BA_t2","BA_pred_calib_lo","BA_pred_calib_hi"]].dropna()
print("\nmean BA interval width as %% of obs:", round(100*((m.BA_pred_calib_hi-m.BA_pred_calib_lo)/m.BA_t2).mean(),1))
