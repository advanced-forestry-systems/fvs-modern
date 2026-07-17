import os,sys,importlib.util,sqlite3
import pandas as pd
CAL="/users/PUOM0008/crsfaaron/fvs-modern/calibration"
sys.path.insert(0,CAL)
spec=importlib.util.spec_from_file_location("H",os.path.join(CAL,"run_5model_scorecard_cohort_vol.py"))
H=importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
BIN="/fs/scratch/PUOM0008/crsfaaron/fvs-modern-volume/bin-r9vol-trace/FVSne"
import subprocess,tempfile
# a few large trees: red spruce(97) and red maple(316) at DBH 16,20,26
trees=[(97,16,60),(97,20,68),(97,26,78),(316,16,55),(316,20,62),(316,26,72)]
for CT in ["F","I"]:
    d=f"/fs/scratch/PUOM0008/crsfaaron/fvs-modern-volume/rescore/trace_{CT}"
    os.makedirs(d,exist_ok=True)
    tf=os.path.join(d,"nvbtrace.txt")
    if os.path.exists(tf): os.remove(tf)
    rows=[]
    for k,(sp,dbh,ht) in enumerate(trees):
        rows.append(dict(stand_id="TR",plot_id=1,tree_id=k+1,tree_count=1.0,species=sp,diameter=float(dbh),ht=float(ht),crratio=50))
    tree_df=pd.DataFrame(rows); stand_df=H.build_standinit("TR",2015,"ne",None)
    db=os.path.join(d,"FVS_Data.db"); con=sqlite3.connect(db)
    stand_df.to_sql("fvs_standinit",con,if_exists="replace",index=False)
    tree_df.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
    key=os.path.join(d,"run.key")
    open(key,"w").write(H.KEYFILE.format(sid="TR",db=db,ncyc=1,calib_kw=""))
    env=os.environ.copy(); env["NVBTRACE"]="1"; env["FVSBFCTYPE"]=CT; env["FVSCFCTYPE"]=CT
    subprocess.run([BIN,f"--keywordfile={key}"],cwd=d,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=120)
    print(f"\n===== TRACE CTYPE={CT} (cols: DBH HTTOT MTOPP HT1PRD STUMP VOL10bf VOL4cf NUMSEG) =====")
    if os.path.exists(tf):
        print(open(tf).read().strip())
    else:
        print("no trace file")
