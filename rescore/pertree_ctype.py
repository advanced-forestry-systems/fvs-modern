import os,sys,importlib.util
import pandas as pd, numpy as np
CAL="/users/PUOM0008/crsfaaron/fvs-modern/calibration"
sys.path.insert(0,CAL)
spec=importlib.util.spec_from_file_location("H",os.path.join(CAL,"run_5model_scorecard_cohort_vol.py"))
H=importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
BIN=sys.argv[1]; CT=sys.argv[2]  # 'F' or 'I'
H.VARIANT_BINS["ne"]=BIN
os.environ["FVSBFCTYPE"]=CT; os.environ["FVSCFCTYPE"]=CT
FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
SWSET={12,43,68,70,71,90,91,94,95,96,97,105,123,125,126,129,130,241,242,299}
cols=["CN","PLT_CN","STATUSCD","SPCD","DIA","HT","CR","TPA_UNADJ","VOLCFNET","VOLBFNET"]
t=pd.read_csv(os.path.join(FIA,"ME_TREE.csv"),usecols=lambda c:c in cols,low_memory=False)
t=t[(t.STATUSCD==1)&(t.DIA>=10)&(t.HT>0)&(t.VOLBFNET>0)&(t.VOLCFNET>0)&(t.TPA_UNADJ>0)].copy()
t["SW"]=t.SPCD.isin(SWSET)
samp=[]
for sw in [True,False]:
    sub=t[t.SW==sw]
    for lo,hi in [(10,13),(13,16),(16,19),(19,23),(23,40)]:
        b=sub[(sub.DIA>=lo)&(sub.DIA<hi)]
        if len(b)>0: samp.append(b.sample(min(20,len(b)),random_state=int(lo)+int(sw)))
S=pd.concat(samp).reset_index(drop=True)
recs=[]
for i,r in S.iterrows():
    sid=f"NE_{i:06d}"; tdf=pd.DataFrame([r]); tdf["PLT_CN"]=r.PLT_CN
    tree_df=H.build_treeinit(tdf,sid)
    if len(tree_df)==0: continue
    stand_df=H.build_standinit(sid,2015,"ne",None)
    res,err=H.run_stand(sid,stand_df,tree_df,"ne","fvs_base",1)
    if err or res is None or len(res)==0: continue
    c0=res.iloc[0]; tpa=float(tree_df.iloc[0]["tree_count"])
    recs.append(dict(SW=bool(r.SW),DBH=float(r.DIA),
        FVS_BdFt=float(c0.get("BdFt",np.nan))/tpa,FIA_VOLBFNET=float(r.VOLBFNET),
        FVS_MCuFt=float(c0.get("MCuFt",np.nan))/tpa,FIA_VOLCFNET=float(r.VOLCFNET)))
D=pd.DataFrame(recs)
D["BF_ratio"]=D.FVS_BdFt/D.FIA_VOLBFNET; D["CF_ratio"]=D.FVS_MCuFt/D.FIA_VOLCFNET
print(f"\n===== CTYPE={CT}  bin={os.path.basename(os.path.dirname(BIN))} =====")
for sw in [True,False]:
    g=D[D.SW==sw]; lab="SW" if sw else "HW"
    print(f"--- {lab} (n={len(g)}) ---  {'DBHbin':<8}{'n':>4}{'BF_ratio':>10}{'CF_ratio':>10}")
    for lo,hi in [(10,13),(13,16),(16,19),(19,23),(23,40)]:
        b=g[(g.DBH>=lo)&(g.DBH<hi)]
        if len(b)==0: continue
        print(f"{'':<19}{str(lo)+'-'+str(hi):<8}{len(b):>4}{b.BF_ratio.mean():>10.3f}{b.CF_ratio.mean():>10.3f}")
# aggregate weighted-ish overall mean ratio for large (>=16)
for sw,lab in [(True,'SW'),(False,'HW')]:
    g=D[(D.SW==sw)&(D.DBH>=16)]
    print(f"OVERALL {lab} DBH>=16: BF_ratio_mean={g.BF_ratio.mean():.3f} CF_ratio_mean={g.CF_ratio.mean():.3f} n={len(g)}")
print("DONE_CT",CT,flush=True)
