import os,sys,math
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P)
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
look={k:v for k,v in zip(brms.key,brms["SDImax.median"])}
print("brms rows",len(brms),"unique keys",len(look))
# NE: join to FIA plots, compare brms(English) to observed SDI
G.VARIANT_STATES["ne"]=(9,23,25,33,36,50); INV={9:"CT",23:"ME",25:"MA",33:"NH",36:"NY",50:"VT"}
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
states=["CT","ME","MA","NH","NY","VT"]
pl=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","STATECD","UNITCD","COUNTYCD","PLOT","MEASYEAR"),low_memory=False) for s in states]
plot=pd.concat(pl,ignore_index=True)
plot["key"]=plot.STATECD.astype(str)+"-"+plot.UNITCD.astype(str)+"-"+plot.COUNTYCD.astype(str)+"-"+plot.PLOT.astype(str)
plot["brms_metric"]=plot.key.map(look)
matched=plot.dropna(subset=["brms_metric"])
print("NE plots",len(plot),"matched to brms",len(matched),"(%.0f%%)"%(100*len(matched)/len(plot)))
# observed SDI for a sample, compare to brms/2.471 (English)
samp=matched.drop_duplicates("CN").sample(n=min(400,len(matched)),random_state=3)
tr=G.load_fia_trees("ne",samp.CN.astype("int64").tolist(),Path(FIA))
rows=[]
for cn,grp in tr.groupby("PLT_CN"):
    d=grp[(grp.STATUSCD==1)&(grp.DIA>0)&(grp.TPA_UNADJ>0)]
    if len(d)<5: continue
    sdi_eng=(d.TPA_UNADJ*(d.DIA/10.0)**1.605).sum()
    rows.append((int(cn),sdi_eng))
obs=pd.DataFrame(rows,columns=["CN","obs_sdi_eng"]); obs["CN"]=obs.CN.astype("int64")
samp["CN"]=samp.CN.astype("int64"); j=samp.merge(obs,on="CN")
j["brms_eng"]=j.brms_metric/2.471
print("NE matched plots with trees:",len(j))
print("observed SDI English: median %.0f p95 %.0f"%(j.obs_sdi_eng.median(),j.obs_sdi_eng.quantile(.95)))
print("brms SDImax English (/2.471): median %.0f p95 %.0f"%(j.brms_eng.median(),j.brms_eng.quantile(.95)))
print("brms metric: median %.0f"%(j.brms_metric.median()))
print("frac plots where obs_sdi < brms_eng (should be ~most):",round((j.obs_sdi_eng<j.brms_eng).mean(),2))
print("DONE_VERIFY")
