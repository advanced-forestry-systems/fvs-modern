#!/usr/bin/env python3
"""ACD (Acadian) and ADK (Adirondacks) as NE sub-variants (customR), calibration layer (2026-06-17).
Run the NE engine on the ACD region (ME,NH,VT) and ADK region (NY Adirondack counties); default vs
calibrated (brms SDImax + ingrowth + BAIMULT). The customR ACD/ADK growth equations are Ben Rice's; this
applies the FIA-grounded calibration layer (SDImax/density/size) that sits above them. Env NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108; SDIc=2.471; VOLc=0.0699055
NSAMP=int(os.environ.get("NSAMP","400")); SEEDR=int(os.environ.get("SEED","5")); FRAC=0.7; BAI=0.90
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/adk_acd.csv"))
FIPS={"ME":23,"NH":33,"VT":50,"NY":36}; INV={v:k for k,v in FIPS.items()}
ENGINE="ne"; G.VARIANT_STATES["ne"]=(23,33,50,36)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
ADK_CO={19,31,33,35,41,43,49,89,113,115}  # NY Adirondack-region county FIPS (Clinton,Essex,Franklin,Fulton,Hamilton,Herkimer,Lewis,StLawrence,Warren,Washington)
INGR_acd=36.1/100.0; INGR_adk=33.0/100.0
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def seed_rows(tdf,rt,spp):
    if rt<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1; per=rt/len(spp); rows=[{"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40} for j,sp in enumerate(spp)]
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
def load_meta(states):
    key={}; cls={}
    for s in states:
        pf=Path(FIA)/f"{s}_PLOT.csv"
        if pf.exists():
            pp=pd.read_csv(pf,usecols=lambda c:c in ("CN","STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False)
            for _,r in pp.iterrows():
                try: key[int(r.CN)]="%d-%d-%d-%d"%(r.STATECD,r.UNITCD,r.COUNTYCD,r.PLOT)
                except: pass
        cf=Path(FIA)/f"{s}_COND.csv"
        if cf.exists():
            c=pd.read_csv(cf,usecols=lambda x:x in ("PLT_CN","TRTCD1","TRTCD2","TRTCD3","DSTRBCD1","DSTRBCD2","DSTRBCD3"),low_memory=False)
            for cn,g2 in c.groupby("PLT_CN"):
                try: cn=int(cn)
                except: continue
                trt=g2[["TRTCD1","TRTCD2","TRTCD3"]].fillna(0).values; dst=g2[["DSTRBCD1","DSTRBCD2","DSTRBCD3"]].fillna(0).values
                cls[cn]="harvested" if (trt==10).any() else ("disturbed" if (dst>0).any() else "undisturbed")
    return key,cls
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def run(std,tdf,sid,kw,yrs):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+ENGINE),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def bias(f,o):
    m=[(x,y) for x,y in zip(f,o) if y>0]; return float("nan") if not m else 100*sum(x-y for x,y in m)/sum(y for _,y in m)
# load NY + ME/NH/VT plots
allstates=["ME","NH","VT","NY"]
pl=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD","COUNTYCD"),low_memory=False) for s in allstates]
plot=pd.concat(pl,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
key,cls=load_meta(allstates)
def region_rows(name):
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)]
    if name=="acd": rem=rem[rem.STATECD.isin([23,33,50])]
    else: rem=rem[(rem.STATECD==36)&(rem.COUNTYCD.isin(ADK_CO))]
    return rem.sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
rows=[]
for name,rate in [("acd",INGR_acd),("adk",INGR_adk)]:
    rem=region_rows(name)
    tr1=G.load_fia_trees(ENGINE,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(ENGINE,rem.CN.tolist(),Path(FIA))
    A={k:{m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")} for k in ("default","calibrated")}; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA,oTPH,oQMD=metr(bl)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,ENGINE); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        sdival=BRMS.get(key.get(cn,""),None); kwc=""
        if sdival and sdival>0: kwc+=sdimax_kw(sdival)+"\n"
        kwc+=allsp("BAIMULT",BAI)+"\n"
        tdf_c=seed_rows(tdf,rate*(yrs/10.0)*float(al.TPA_UNADJ.sum())*FRAC,list(al.SPCD.value_counts().index[:3]))
        sd=run(std,tdf,sid,"",yrs); sc=run(std,tdf_c,sid,kwc,yrs)
        if sd is None or len(sd)==0 or sc is None or len(sc)==0: continue
        nu+=1
        for k,s in [("default",sd),("calibrated",sc)]:
            l=s.iloc[-1]
            A[k]["BA"]["f"].append(g(l,"BA")*M2HA); A[k]["BA"]["o"].append(oBA)
            A[k]["TPH"]["f"].append(g(l,"Tpa")*TPHc); A[k]["TPH"]["o"].append(oTPH)
            A[k]["QMD"]["f"].append(g(l,"QMD")*CMc); A[k]["QMD"]["o"].append(oQMD)
    if nu<8: print(name,"too few (n=%d)"%nu); sys.stdout.flush(); continue
    row={"subvariant":name,"engine":ENGINE,"n":nu,
         "def_BA%":round(bias(A["default"]["BA"]["f"],A["default"]["BA"]["o"]),1),"cal_BA%":round(bias(A["calibrated"]["BA"]["f"],A["calibrated"]["BA"]["o"]),1),
         "def_TPH%":round(bias(A["default"]["TPH"]["f"],A["default"]["TPH"]["o"]),1),"cal_TPH%":round(bias(A["calibrated"]["TPH"]["f"],A["calibrated"]["TPH"]["o"]),1),
         "def_QMD%":round(bias(A["default"]["QMD"]["f"],A["default"]["QMD"]["o"]),1),"cal_QMD%":round(bias(A["calibrated"]["QMD"]["f"],A["calibrated"]["QMD"]["o"]),1)}
    rows.append(row)
    print("%-4s (NE+customR) n=%-3d | BA %+5.1f->%+5.1f | TPH %+6.1f->%+6.1f | QMD %+5.1f->%+5.1f"%(name,nu,row["def_BA%"],row["cal_BA%"],row["def_TPH%"],row["cal_TPH%"],row["def_QMD%"],row["cal_QMD%"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_ADK")
