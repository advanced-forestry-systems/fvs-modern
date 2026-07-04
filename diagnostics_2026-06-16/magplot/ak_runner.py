#!/usr/bin/env python3
"""Run FVS-AK (default) on coastal-BC MAGPlot stands and compare to observed t2.
Adapted from magplot_fvs_runner.py: standalone binary + SQLite keyfile, one stand/subprocess."""
import os, sqlite3, subprocess, tempfile, argparse, math, sys
import pandas as pd, numpy as np
MP="/users/PUOM0008/crsfaaron/magplot"; ROOT="/users/PUOM0008/crsfaaron/fvs-modern"
LIB=os.path.join(ROOT,"lib")
ACRES_PER_HA=2.4710538147; INCH_PER_CM=1/2.54; FT_PER_M=1/0.3048
# BC genus.epithet -> FIA SPCD (coastal + common interior)
GS={"PINU.CON":108,"PSEU.MEN":202,"TSUG.HET":263,"THUJ.PLI":242,"PICE.GLA":94,
"ABIE.LAS":19,"POPU.TRE":746,"LARI.OCC":73,"BETU.PAP":375,"ABIE.AMA":11,
"PICE.SIT":98,"PICE.ENG":93,"ALNU.RUB":351,"PICE.MAR":95,"SALI.SPP":920,
"PINU.MON":119,"PINU.PON":122,"POPU.BAL":741,"TSUG.MER":264,"CHAM.NOO":42,
"CALL.NOO":42,"ABIE.GRA":17,"POPU.TRI":747,"ACER.MAC":312,"ARBU.MEN":361,
"QUER.GAR":815,"PINU.ALB":101,"TAXU.BRE":231,"CORN.NUT":492,"PRUN.SPP":760,
"BETU.SPP":370,"ACER.SPP":310,"POPU.SPP":740,"PICE.SPP":90,"ABIE.SPP":10,
"PINU.SPP":100,"TSUG.SPP":260}
KEY="""STDIDENT
{sid}
DATABASE
DSNIN
{db}
DSNOUT
{db}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'
ENDSQL
END
DATABASE
SUMMARY            2
END
TIMEINT            0         {clen}
NUMCYCLE          {ncyc}
PROCESS
STOP
"""
def standinit(sid,yr,lat,lon,elevft):
    return pd.DataFrame([{"stand_id":sid,"variant":"AK","inv_year":int(yr),
    "latitude":round(float(lat),2),"longitude":round(float(lon),2),"region":10,"forest":5,"district":0,
    "basal_area_factor":0.0,"inv_plot_size":1.0,"brk_dbh":5.0,"num_plots":1,"age":80,
    "aspect":0,"slope":10,"elevft":int(elevft),"site_species":98,"site_index":80,
    "state":2,"county":1,"forest_type":301,"sam_wt":1.0}])
def treeinit(g,sid):
    rows=[]; tid=1
    for _,r in g.iterrows():
        sp=GS.get(str(r["species_gs"]).strip().upper())
        if sp is None: continue
        try: d=float(r["dbh"]); s=float(r["stem_ha"])
        except: continue
        if not(np.isfinite(d) and d>0 and np.isfinite(s) and s>0): continue
        ht=0.0
        try:
            h=float(r["height"]);
            if np.isfinite(h) and h>0: ht=h*FT_PER_M
        except: pass
        rows.append({"stand_id":sid,"plot_id":1,"tree_id":tid,
          "tree_count":round(s/ACRES_PER_HA,5),"species":sp,
          "diameter":round(d*INCH_PER_CM,3),"ht":round(ht,1),"crratio":40}); tid+=1
    return pd.DataFrame(rows)
def run_stand(sid,tdf,yr,lat,lon,elevft,interval):
    binary=os.path.join(LIB,"FVSak")
    ncyc=max(1,math.ceil(interval/10)); clen=max(1,round(interval/ncyc))
    with tempfile.TemporaryDirectory() as d:
        db=os.path.join(d,"FVS_Data.db"); con=sqlite3.connect(db)
        standinit(sid,yr,lat,lon,elevft).to_sql("fvs_standinit",con,if_exists="replace",index=False)
        tdf.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
        key=os.path.join(d,"ak.key"); open(key,"w").write(KEY.format(sid=sid,db=db,clen=clen,ncyc=ncyc))
        try:
            subprocess.run([binary,f"--keywordfile={key}"],cwd=d,stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,timeout=120)
            con=sqlite3.connect(db); df=pd.read_sql_query("SELECT * FROM FVS_Summary2",con); con.close()
            return df,clen*ncyc,None
        except Exception as e: return None,0,str(e)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--nstands",type=int,default=30)
    ap.add_argument("--min-trees",type=int,default=8); ap.add_argument("--coastal",type=float,default=-123.0)
    ap.add_argument("--out",default=f"{MP}/ak_fvs_results.csv"); a=ap.parse_args()
    sites=pd.read_csv(f"{MP}/magp_sites.csv",low_memory=False); sites["magp_site_id"]=sites["magp_site_id"].astype(str)
    coastal=set(sites.loc[(sites["province"]=="BC")&(sites["longitude"]<=a.coastal),"magp_site_id"])
    print(f"coastal BC sites (lon<={a.coastal}): {len(coastal)}",flush=True)
    pairs=pd.read_csv(f"{MP}/magplot_ak_bc_pairs_clean.csv"); pairs["SITE"]=pairs["SITE"].astype(str)
    pairs=pairs[pairs["SITE"].isin(coastal)].copy()
    print(f"coastal pairs: {len(pairs)}",flush=True)
    tr=pd.read_csv(f"{MP}/magp_trees_bc.csv",low_memory=False); tr["magp_site_id"]=tr["magp_site_id"].astype(str)
    tr["species_gs"]=tr["species_g"].astype(str)+"."+tr["species_s"].astype(str)
    tr=tr[(tr["tree_status"].astype(str).str.startswith("L"))&(tr["dbh"]>0)&(tr["stem_ha"]>0)]
    sll=sites.set_index("magp_site_id")[["latitude","longitude","elevation"]].to_dict("index")
    out=[]; done=0; cov_acc=[]
    for _,p in pairs.iterrows():
        site=p["SITE"]; ptype="site"
        g=tr[(tr["magp_site_id"]==site)&(tr["meas_num"]==int(p["meas1"]))]
        if not len(g): continue
        if len(g)<a.min_trees: continue
        cov=g["species_gs"].str.upper().isin(GS).mean(); cov_acc.append(cov)
        tdf=treeinit(g,site)
        if len(tdf)<a.min_trees: continue
        meta=sll.get(site,{"latitude":52,"longitude":-126,"elevation":300})
        elevft=int((meta.get("elevation") or 300)*FT_PER_M)
        df,projyrs,err=run_stand(site,tdf,int(p["y1"]),meta["latitude"],meta["longitude"],elevft,int(p["interval_years"]))
        if err or df is None or not len(df):
            print(f"  {site} {ptype}: ERR {err}",flush=True); continue
        yn=df.iloc[-1]
        ba_obs_m2ha=p["BA_t2_obs"]
        out.append(dict(site=site,ptype=ptype,ntrees=len(tdf),interval=p["interval_years"],projyrs=projyrs,
            BA_proj_t1_m2ha=df.iloc[0]["BA"]/4.356,BA_proj_m2ha=yn["BA"]/4.356,BA_t1_obs_m2ha=p["BA_t1_m2ha"],BA_obs_m2ha=ba_obs_m2ha,Tpa_proj=yn["Tpa"],TPH_obs=p["TPH_t2_obs"],
            QMD_proj=yn.get("QMD",np.nan),QMD_obs_cm=p["QMD_t2_obs"]))
        print(f"  {site} {ptype}: {len(tdf)}tr cov{cov*100:.0f}% proj{projyrs}y BAm2 {df.iloc[0]['BA']/4.356:.1f}->{yn['BA']/4.356:.1f} obs {ba_obs_m2ha:.1f}",flush=True)
        done+=1
        if done>=a.nstands: break
    print(f"mean crosswalk coverage: {100*np.mean(cov_acc):.1f}%",flush=True)
    if out:
        R=pd.DataFrame(out); R.to_csv(a.out,index=False)
        levbias=100*(R["BA_proj_m2ha"]-R["BA_obs_m2ha"]).sum()/R["BA_obs_m2ha"].sum()
        R["dproj"]=R["BA_proj_m2ha"]-R["BA_proj_t1_m2ha"]; R["dobs"]=R["BA_obs_m2ha"]-R["BA_t1_obs_m2ha"]
        incbias=100*(R["dproj"].sum()-R["dobs"].sum())/R["dobs"].sum() if R["dobs"].sum()!=0 else float("nan")
        print(f"AK n={len(R)}: LEVEL BA bias {levbias:+.1f}% (proj {R['BA_proj_m2ha'].mean():.1f} vs obs {R['BA_obs_m2ha'].mean():.1f})",flush=True)
        print(f"AK n={len(R)}: INCREMENT BA bias {incbias:+.1f}% (proj dBA {R['dproj'].mean():.2f} vs obs dBA {R['dobs'].mean():.2f} m2/ha/{R['interval'].mean():.0f}y)",flush=True)
        print(f"  proj t1 {R['BA_proj_t1_m2ha'].mean():.1f} vs obs t1 {R['BA_t1_obs_m2ha'].mean():.1f} m2/ha (tag-limit offset)",flush=True)
    print("DONE_AK_PILOT",flush=True)
main()
