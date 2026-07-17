import pandas as pd, glob, os
COND="/fs/scratch/PUOM0008/crsfaaron/FIA"
need=["PLT_CN","COND_STATUS_CD","STDORGCD","CONDPROP_UNADJ"]
rows=[]; skipped=[]
for f in sorted(glob.glob(f"{COND}/*_COND.csv")):
    st=os.path.basename(f).split("_")[0]
    c=pd.read_csv(f,nrows=0).columns.tolist()
    if not all(x in c for x in need):
        skipped.append(st); continue
    d=pd.read_csv(f,usecols=need,dtype={"PLT_CN":str},low_memory=False)
    d=d[d.COND_STATUS_CD==1]
    pl=d.assign(pl_prop=(d.STDORGCD==1)*d.CONDPROP_UNADJ).groupby("PLT_CN").agg(
        pl=("pl_prop","sum"),tot=("CONDPROP_UNADJ","sum"))
    pl["plantation"]=(pl.pl>0.5*pl.tot).astype(int)
    pf=pl.reset_index()[["PLT_CN","plantation"]]; pf["STATE"]=st
    rows.append(pf)
allp=pd.concat(rows,ignore_index=True)
allp["PLT_CN"]=allp.PLT_CN.str.replace(r"\.0$","",regex=True)
allp.to_csv("plt_plantation.csv",index=False)
print("plots",len(allp),"plantation",int(allp.plantation.sum()),
      f"{100*allp.plantation.mean():.1f}%  skipped:",skipped)
