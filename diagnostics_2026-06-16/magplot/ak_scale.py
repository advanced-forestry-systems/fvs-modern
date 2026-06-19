import os, sys
import pandas as pd, numpy as np, math
sys.path.insert(0,"/users/PUOM0008/crsfaaron/magplot")
os.environ["FVS_LIB_DIR"]="/users/PUOM0008/crsfaaron/fvs-modern/lib"
import magplot_fvs_runner as R
MP="/users/PUOM0008/crsfaaron/magplot"
R.GS_TO_SPCD={**R.GS_TO_SPCD,"TSUG.HET":263,"THUJ.PLI":242,"PICE.SIT":98,"PSEU.MEN":202,"ABIE.AMA":11,"PINU.CON":108,"PICE.GLA":94,"ABIE.LAS":19,"PICE.ENG":93,"PICE.MAR":95,"ALNU.RUB":351,"POPU.TRE":746,"LARI.OCC":73,"BETU.PAP":375}
bc=pd.read_csv(f"{MP}/magp_trees_bc.csv",low_memory=False); bc["magp_site_id"]=bc["magp_site_id"].astype(str)
bc["species_gs"]=bc["species_g"].astype(str)+"."+bc["species_s"].astype(str)
g=bc[(bc["magp_site_id"]=="30353133")&(bc["tree_status"].astype(str).str.startswith("L"))&(bc["dbh"]>0)&(bc["stem_ha"]>0)]
g=g[g["meas_num"]==g["meas_num"].max()]
tdf=R.build_treeinit_magplot(g,"BCX")
print("base sum tree_count",round(tdf['tree_count'].sum(),1),"n",len(tdf))
for scale in [1,10]:
    t2=tdf.copy(); t2["tree_count"]=t2["tree_count"]*scale
    df,err=R.run_stand("BCX",t2,"ak","default",1)
    y=df.iloc[0]; print(f"  scale x{scale}: FVS TPA {y['Tpa']:.1f}  BA {y['BA']:.1f}ft2")
# also: collapse to ONE record per species with summed tree_count (avoid many small records)
agg=g.copy()
agg["sp"]=agg["species_gs"].map(R.GS_TO_SPCD)
rows=[]; tid=1
for sp,sub in agg.groupby("sp"):
    # weighted mean dbh, total count
    w=sub["stem_ha"]; mdbh=np.sqrt(np.sum(w*sub["dbh"]**2)/np.sum(w))
    rows.append({"stand_id":"BCX","plot_id":1,"tree_id":tid,"tree_count":round(sub["stem_ha"].sum()/2.471,4),
      "species":int(sp),"diameter":round(mdbh/2.54,3),"ht":0.0,"crratio":40}); tid+=1
import pandas as pd
t3=pd.DataFrame(rows); print("collapsed records:",len(t3),"sum TPA",round(t3['tree_count'].sum(),1))
df,err=R.run_stand("BCX",t3,"ak","default",1)
if not err: y=df.iloc[0]; print(f"  collapsed: FVS TPA {y['Tpa']:.1f} BA {y['BA']:.1f}ft2 ={y['BA']/4.356:.1f}m2/ha")
