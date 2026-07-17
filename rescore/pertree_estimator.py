import os,sys,importlib.util,math,sqlite3,tempfile,subprocess,shutil
import pandas as pd, numpy as np
CAL="/users/PUOM0008/crsfaaron/fvs-modern/calibration"
sys.path.insert(0,CAL)
spec=importlib.util.spec_from_file_location("H",os.path.join(CAL,"run_5model_scorecard_cohort_vol.py"))
H=importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
BIN="/fs/scratch/PUOM0008/crsfaaron/fvs-modern-volume/bin-r9vol/FVSne"
H.VARIANT_BINS["ne"]=BIN
FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
SWSET={12,43,68,70,71,90,91,94,95,96,97,105,123,125,126,129,130,241,242,299} # softwood SPCD (spruce,fir,pine,cedar,hemlock,larch,etc)
# load ME live large trees
cols=["CN","PLT_CN","STATUSCD","SPCD","DIA","HT","CR","TPA_UNADJ","VOLCFNET","VOLBFNET"]
t=pd.read_csv(os.path.join(FIA,"ME_TREE.csv"),usecols=lambda c:c in cols,low_memory=False)
t=t[(t.STATUSCD==1)&(t.DIA>=10)&(t.HT>0)&(t.VOLBFNET>0)&(t.VOLCFNET>0)&(t.TPA_UNADJ>0)].copy()
t["SW"]=t.SPCD.isin(SWSET)
# sample across DBH bins, both SW and HW
rng=np.random.default_rng(42)
samp=[]
for sw in [True,False]:
    sub=t[t.SW==sw]
    for lo,hi in [(10,13),(13,16),(16,19),(19,23),(23,40)]:
        b=sub[(sub.DIA>=lo)&(sub.DIA<hi)]
        if len(b)>0:
            samp.append(b.sample(min(20,len(b)),random_state=int(lo)+int(sw)))
S=pd.concat(samp).reset_index(drop=True)
print("sampled",len(S),"trees",flush=True)
recs=[]
for i,r in S.iterrows():
    sid=f"NE_{i:06d}"
    tdf=pd.DataFrame([r]); tdf["PLT_CN"]=r.PLT_CN
    tree_df=H.build_treeinit(tdf,sid)
    if len(tree_df)==0: continue
    stand_df=H.build_standinit(sid,2015,"ne",None)
    res,err=H.run_stand(sid,stand_df,tree_df,"ne","fvs_base",1)
    if err or res is None or len(res)==0:
        continue
    c0=res.iloc[0]
    tpa=float(tree_df.iloc[0]["tree_count"])
    bdft_tree=float(c0.get("BdFt",np.nan))/tpa
    mcuft_tree=float(c0.get("MCuFt",np.nan))/tpa
    recs.append(dict(CN=r.CN,SPCD=int(r.SPCD),SW=bool(r.SW),DBH=float(r.DIA),HT=float(r.HT),
        TPA=tpa,FVS_BdFt=bdft_tree,FIA_VOLBFNET=float(r.VOLBFNET),
        FVS_MCuFt=mcuft_tree,FIA_VOLCFNET=float(r.VOLCFNET)))
D=pd.DataFrame(recs)
D["BF_ratio"]=D.FVS_BdFt/D.FIA_VOLBFNET
D["CF_ratio"]=D.FVS_MCuFt/D.FIA_VOLCFNET
D["BF_per_CF_FVS"]=D.FVS_BdFt/D.FVS_MCuFt
D["BF_per_CF_FIA"]=D.FIA_VOLBFNET/D.FIA_VOLCFNET
D.to_csv("/fs/scratch/PUOM0008/crsfaaron/fvs-modern-volume/rescore/pertree_estimator.csv",index=False)
print("\n=== Per-tree estimator (cycle-0, identical inputs) FVS NSVB vs FIA NSVB ===",flush=True)
for sw in [True,False]:
    g=D[D.SW==sw]
    lab="SW" if sw else "HW"
    print(f"\n--- {lab} (n={len(g)}) ---")
    print(f"{'DBHbin':<10}{'n':>4}{'BF_ratio':>10}{'CF_ratio':>10}{'BFpCF_FVS':>11}{'BFpCF_FIA':>11}")
    for lo,hi in [(10,13),(13,16),(16,19),(19,23),(23,40)]:
        b=g[(g.DBH>=lo)&(g.DBH<hi)]
        if len(b)==0: continue
        print(f"{str(lo)+'-'+str(hi):<10}{len(b):>4}{b.BF_ratio.mean():>10.3f}{b.CF_ratio.mean():>10.3f}{b.BF_per_CF_FVS.mean():>11.2f}{b.BF_per_CF_FIA.mean():>11.2f}")
print("\nDONE_PERTREE",flush=True)
