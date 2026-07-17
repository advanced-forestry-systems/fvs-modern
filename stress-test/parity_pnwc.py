#!/usr/bin/env python3
import os,sys,numpy as np,pandas as pd
PR="/users/PUOM0008/crsfaaron/fvs-modern"
sys.path.insert(0,PR); sys.path.insert(0,"/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
TREEDIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h"
EQDIR=os.path.join(SCR,"conus_eq_proj","out_conus_eq")
BA_CONST=np.pi/(4*144); FIPS={41:"OR",53:"WA"}
def fix_si(sdf):
    sdf["inv_plot_size"]=1.0; sdf["num_plots"]=1; sdf["brk_dbh"]=999.0; sdf["basal_area_factor"]=0.0; return sdf
def run_variant(v,ncap=60):
    si=pd.read_csv(os.path.join(SCR,"standinit_by_variant",f"standinit_{v}.csv"),low_memory=False)
    si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
    eqf=os.path.join(EQDIR,f"conus_eq_{v.lower()}_conus_b2_metrics.csv")
    eq=pd.read_csv(eqf); eq=eq[eq["PROJ_YEAR"]==0]; eq["STAND_CN"]=eq["STAND_CN"].astype(str)
    eqm=dict(zip(eq["STAND_CN"],eq["BA_FT2AC"]))
    tcache={}
    res=[]
    done=0
    for _,stand in si.iterrows():
        if done>=ncap: break
        cn=stand["STAND_CN"]
        try: state=FIPS[int(float(stand["STATE"]))]
        except: continue
        if state not in tcache:
            tt=pd.read_csv(os.path.join(TREEDIR,f"{state}_FVS_TREEINIT_PLOT.csv"),low_memory=False)
            tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
            tcache[state]={k:vv for k,vv in tt.groupby("STAND_CN")}
        fr_rows=tcache[state].get(cn)
        if fr_rows is None or fr_rows.empty: continue
        # observed
        fr2=fr_rows[(fr_rows["DIAMETER"]>=1.0)&(fr_rows["TREE_COUNT"]>0)]
        if fr2.empty: continue
        obs=(BA_CONST*fr2["DIAMETER"]**2*fr2["TREE_COUNT"]).sum()
        sid=f"S{cn}"
        plot_data={"INVYR":int(float(stand.get("INV_YEAR") or 2010)),"LAT":stand.get("LATITUDE"),"LON":stand.get("LONGITUDE"),
                   "ELEV":stand.get("ELEVFT") or 500,"SLOPE":stand.get("SLOPE") or 10,"ASPECT":stand.get("ASPECT") or 180,"STDAGE":stand.get("AGE") or 50}
        sdf=fix_si(P.build_fvs_standinit(plot_data,sid,v.lower()))
        recs=[]
        for i,t in enumerate(fr_rows.itertuples(index=False)):
            d=getattr(t,"DIAMETER",np.nan)
            if pd.isna(d) or float(d)<1.0: continue
            recs.append({"stand_id":sid,"plot_id":1,"tree_id":i+1,"tree_count":float(getattr(t,"TREE_COUNT",1.0) or 1.0),
                         "species":int(float(getattr(t,"SPECIES",0) or 0)),"diameter":round(float(d),1),
                         "ht":(round(float(getattr(t,"HT",0)),0) if pd.notna(getattr(t,"HT",0)) and float(getattr(t,"HT",0) or 0)>0 else 0),"crratio":(int(float(getattr(t,"CRRATIO",0))) if pd.notna(getattr(t,"CRRATIO",0)) and float(getattr(t,"CRRATIO",0) or 0)>0 else 0)})
        tdf=pd.DataFrame(recs)
        if tdf.empty: continue
        try:
            fr=P.run_fvs_projection(sdf,tdf,sid,v.lower(),config_version=None,num_cycles=1,cycle_length=5)
        except Exception as e: continue
        tls=fr["treelists"]
        if not tls: res.append((cn,obs,0.0,eqm.get(cn,np.nan))); done+=1; continue
        tl=tls[min(tls.keys())]
        tc="TPA" if "TPA" in tl.columns else "Tpa"; dc="DBH" if "DBH" in tl.columns else "Dbh"
        eng=(BA_CONST*tl[dc].astype(float)**2*tl[tc].astype(float)).sum()
        res.append((cn,obs,eng,eqm.get(cn,np.nan))); done+=1
    df=pd.DataFrame(res,columns=["cn","obs","eng","eq"])
    df=df[df["obs"]>0]
    r=df["eng"]/df["obs"]; within5=np.mean(np.abs(df["eng"]-df["obs"])/df["obs"]<=0.05)*100
    ba0=np.mean(df["eng"]==0)*100
    eng_eq=np.nanmean(df["eng"]/df["eq"].replace(0,np.nan))
    print(f"{v}: n={len(df)} obsBA={df['obs'].mean():.2f} engV3BA={df['eng'].mean():.2f} eqBA={df['eq'].mean():.2f} eng/obs={r.mean():.3f} eng/eq={eng_eq:.3f} within5%={within5:.1f}% BA0%={ba0:.1f}%")
for v in ["PN","WC"]: run_variant(v,60)
