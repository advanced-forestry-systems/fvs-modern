import os,sys,sqlite3,subprocess,tempfile,math
import pandas as pd,numpy as np
from multiprocessing import Pool
sys.path.insert(0,"/users/PUOM0008/crsfaaron/magplot")
LIB="/users/PUOM0008/crsfaaron/fvs-modern/lib"; MP="/users/PUOM0008/crsfaaron/magplot"; W="/fs/scratch/PUOM0008/crsfaaron/akwork"
import magplot_fvs_runner as R
R.GS_TO_SPCD={**R.GS_TO_SPCD,"TSUG.HET":263,"THUJ.PLI":242,"PICE.SIT":98,"PSEU.MEN":202,"ABIE.AMA":11,"PINU.CON":108,"PICE.GLA":94,"ABIE.LAS":19,"PICE.ENG":93,"PICE.MAR":95,"ALNU.RUB":351,"POPU.TRE":746,"LARI.OCC":73,"BETU.PAP":375,"TSUG.MER":264,"CHAM.NOO":42,"ABIE.GRA":17,"PINU.MON":119,"POPU.BAL":741}
Q=chr(39)
KEY=("STDIDENT\n{sid}\nDATABASE\nDSNIN\n{db}\nDSNOUT\n{db}\nSTANDSQL\n"
 "SELECT * FROM fvs_standinit WHERE stand_id = "+Q+"%StandID%"+Q+"\nENDSQL\nTREESQL\n"
 "SELECT * FROM fvs_treeinit WHERE stand_id = "+Q+"%StandID%"+Q+"\nENDSQL\nEND\n"
 "DATABASE\nSUMMARY            2\nEND\n{bai}TIMEINT            0        {clen}\nNUMCYCLE          {ncyc}\nPROCESS\nSTOP\n")
MULTS=[1.0,2.0,3.0,4.0,6.0]
def run(tdf,sid,ps,ncyc,clen,m):
    d=tempfile.mkdtemp(dir="/tmp"); db=d+"/x.db"; con=sqlite3.connect(db)
    si=R.build_standinit(sid,1990,"ak"); si["inv_plot_size"]=ps
    si.to_sql("fvs_standinit",con,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
    bai="" if m==1.0 else "BAIMULT           0         0{:>10.3f}\n".format(m)
    open(d+"/x.key","w").write(KEY.format(sid=sid,db=db,clen=clen,ncyc=ncyc,bai=bai))
    subprocess.run([LIB+"/FVSak","--keywordfile="+d+"/x.key"],cwd=d,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=120)
    con=sqlite3.connect(db); s=pd.read_sql_query("SELECT Tpa,BA FROM FVS_Summary2",con); con.close(); return s
BC=pd.read_csv(MP+"/magp_trees_bc.csv",low_memory=False); BC["magp_site_id"]=BC["magp_site_id"].astype(str)
BC["species_gs"]=BC["species_g"].astype(str)+"."+BC["species_s"].astype(str)
def work(p):
    site=p["SITE"]; g=BC[(BC["magp_site_id"]==site)&(BC["meas_num"]==int(p["meas1"]))&(BC["tree_status"].astype(str).str.startswith("L"))&(BC["dbh"]>0)&(BC["stem_ha"]>0)]
    if len(g)<18: return None
    tdf=R.build_treeinit_magplot(g,site); tgt=tdf["tree_count"].sum()
    iv=int(p["interval_years"]); ncyc=max(1,math.ceil(iv/10)); clen=max(1,round(iv/ncyc))
    s1=run(tdf,site,1.0,ncyc,clen,1.0)
    if s1 is None or s1.iloc[0]["Tpa"]<=0: return None
    ps=tgt/s1.iloc[0]["Tpa"]
    out={"site":site,"ecoL1":p["ecoL1"],"iv":iv,"t1o":p["BA_t1_m2ha"],"t2o":p["BA_t2_obs"]}
    for m in MULTS:
        s=run(tdf,site,ps,ncyc,clen,m)
        if s is None: return None
        if m==1.0:
            t0=s.iloc[0]["BA"]/4.356
            if p["BA_t1_m2ha"]<=0 or abs(t0-p["BA_t1_m2ha"])/p["BA_t1_m2ha"]>0.30: return None
            out["t0"]=t0
        out["m%.0f"%m]=s.iloc[-1]["BA"]/4.356
    return out
if __name__=="__main__":
    pairs=pd.read_csv(MP+"/magplot_ak_bc_pairs_clean.csv"); pairs["SITE"]=pairs["SITE"].astype(str)
    eco=pd.read_csv(W+"/bc_site_ecoregion.csv"); eco["magp_site_id"]=eco["magp_site_id"].astype(str)
    pairs=pairs.merge(eco,left_on="SITE",right_on="magp_site_id",how="left")
    pairs=pairs[pairs["ecoL1"].notna()].sample(frac=1,random_state=7)
    cand=pd.concat([d.head(400) for _,d in pairs.groupby("ecoL1")])
    recs=[r for r in Pool(4).map(work,[r for _,r in cand.iterrows()]) if r]
    R2=pd.DataFrame(recs); R2.to_csv(W+"/ak_calib_results.csv",index=False)
    print("kept",len(R2),flush=True)
    def inc_bias(df,col):
        dpj=df[col]-df["t0"]; dob=df["t2o"]-df["t1o"]
        return 100*(dpj.sum()-dob.sum())/dob.sum() if dob.sum()!=0 else float("nan")
    print("=== AK increment bias: default vs BAIMULT sweep, by ecoregion ===")
    for e1,d in R2.groupby("ecoL1"):
        if len(d)<8: continue
        biases={m:inc_bias(d,"m%.0f"%m) for m in MULTS}
        best=min(MULTS,key=lambda m:abs(biases[m]))
        print("  %-32s n=%3d  default %+6.1f%%  best BAIMULT %.0fx -> %+6.1f%%"%(str(e1)[:30],len(d),biases[1.0],best,biases[best]),flush=True)
    print("DONE_CALIB")
