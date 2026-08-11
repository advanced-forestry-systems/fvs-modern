import os,sys,sqlite3,subprocess,tempfile,math
import pandas as pd,numpy as np
from multiprocessing import Pool
sys.path.insert(0,"/users/PUOM0008/crsfaaron/magplot")
LIB="/users/PUOM0008/crsfaaron/fvs-modern/lib"; MP="/users/PUOM0008/crsfaaron/magplot"; W="/fs/scratch/PUOM0008/crsfaaron/akwork"
import magplot_fvs_runner as R
R.GS_TO_SPCD={**R.GS_TO_SPCD,"TSUG.HET":263,"THUJ.PLI":242,"PICE.SIT":98,"PSEU.MEN":202,"ABIE.AMA":11,"PINU.CON":108,"PICE.GLA":94,"ABIE.LAS":19,"PICE.ENG":93,"PICE.MAR":95,"ALNU.RUB":351,"POPU.TRE":746,"LARI.OCC":73,"BETU.PAP":375,"TSUG.MER":264,"CHAM.NOO":42,"ABIE.GRA":17,"PINU.MON":119,"POPU.BAL":741}
BC=pd.read_csv(MP+"/magp_trees_bc.csv",low_memory=False); BC["magp_site_id"]=BC["magp_site_id"].astype(str)
BC["species_gs"]=BC["species_g"].astype(str)+"."+BC["species_s"].astype(str)
def once(tdf,sid,ps,ncyc,clen):
    d=tempfile.mkdtemp(dir="/tmp"); db=d+"/x.db"; con=sqlite3.connect(db)
    si=R.build_standinit(sid,1990,"ak"); si["inv_plot_size"]=ps
    si.to_sql("fvs_standinit",con,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
    open(d+"/x.key","w").write(R.KEYFILE.format(sid=sid,db=db,clen=clen,ncyc=ncyc,calib="** DEF"))
    subprocess.run([LIB+"/FVSak","--keywordfile="+d+"/x.key"],cwd=d,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=120)
    con=sqlite3.connect(db)
    try: s=pd.read_sql_query("SELECT Year,Tpa,BA FROM FVS_Summary2",con); con.close(); return s
    except: con.close(); return None
def work(p):
    site=p["SITE"]; g=BC[(BC["magp_site_id"]==site)&(BC["meas_num"]==int(p["meas1"]))&(BC["tree_status"].astype(str).str.startswith("L"))&(BC["dbh"]>0)&(BC["stem_ha"]>0)]
    if len(g)<30: return None
    tdf=R.build_treeinit_magplot(g,site); tgt=tdf["tree_count"].sum()
    iv=int(p["interval_years"]); ncyc=max(1,math.ceil(iv/10)); clen=max(1,round(iv/ncyc))
    s1=once(tdf,site,1.0,ncyc,clen)
    if s1 is None or s1.iloc[0]["Tpa"]<=0: return None
    ps=tgt/s1.iloc[0]["Tpa"]; s=once(tdf,site,ps,ncyc,clen) if abs(ps-1)>0.02 else s1
    if s is None: return None
    t0=s.iloc[0]["BA"]/4.356; tn=s.iloc[-1]["BA"]/4.356
    if p["BA_t1_m2ha"]<=0 or abs(t0-p["BA_t1_m2ha"])/p["BA_t1_m2ha"]>0.15: return None
    return dict(site=site,ecoL1=p["ecoL1"],ecoL3=p["ecoL3"],iv=iv,t0=t0,t1o=p["BA_t1_m2ha"],t2=tn,t2o=p["BA_t2_obs"])
if __name__=="__main__":
    pairs=pd.read_csv(MP+"/magplot_ak_bc_pairs_clean.csv"); pairs["SITE"]=pairs["SITE"].astype(str)
    eco=pd.read_csv(W+"/bc_site_ecoregion.csv"); eco["magp_site_id"]=eco["magp_site_id"].astype(str)
    pairs=pairs.merge(eco,left_on="SITE",right_on="magp_site_id",how="left")
    pairs=pairs[pairs["ecoL1"].notna()].sample(frac=1,random_state=5)
    # balance: up to 250 candidates per major ecoregion
    cand=pd.concat([d.head(250) for _,d in pairs.groupby("ecoL1")])
    recs=[r for r in Pool(4).map(work,[r for _,r in cand.iterrows()]) if r]
    R2=pd.DataFrame(recs); R2.to_csv(W+"/ak_eco_validate_results.csv",index=False)
    print("kept",len(R2),flush=True)
    def rep(df,lbl):
        if len(df)<8: print("  %-34s n=%d (few)"%(lbl,len(df))); return
        lev=100*(df["t2"]-df["t2o"]).sum()/df["t2o"].sum()
        dpj=df["t2"]-df["t0"]; dob=df["t2o"]-df["t1o"]
        inc=100*(dpj.sum()-dob.sum())/dob.sum() if dob.sum()!=0 else float("nan")
        print("  %-32s n=%3d ivl=%2.0fy LEVEL %+5.1f%% INCR %+6.1f%% (projdBA %.2f obs %.2f)"%(lbl,len(df),df["iv"].mean(),lev,inc,dpj.mean(),dob.mean()),flush=True)
    print("=== AK default bias by Level I ecoregion ===")
    for e1,d in R2.groupby("ecoL1"): rep(d,str(e1))
    print("=== Marine West Coast Forest by Level III ===")
    for e3,d in R2[R2.ecoL1=="MARINE WEST COAST FOREST"].groupby("ecoL3"): rep(d,str(e3)[:30])
    print("DONE_ECO_VAL")
