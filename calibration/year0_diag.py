#!/usr/bin/env python3
"""Year-0 volume diagnostic: FVS yr0 vs FIA T1 observed."""
import os, sys, sqlite3, subprocess, tempfile, shutil, math
import pandas as pd
import numpy as np

PROJECT_ROOT = "/users/PUOM0008/crsfaaron/fvs-modern"
LIB_DIR      = os.path.join(PROJECT_ROOT, "src-converted", "bin")
FIA_DIR      = "/fs/scratch/PUOM0008/crsfaaron/FIA"

KEYFILE_TMPL = """\
STDIDENT
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
TIMEINT            0         1
NUMCYCLE          {ncyc}
** DEFAULT PARAMETERS
PROCESS
STOP
"""

STATE_ABBR = {
    23:'ME',24:'MD',25:'MA',33:'NH',34:'NJ',
    36:'NY',42:'PA',50:'VT',54:'WV',9:'CT'
}

def build_standinit(sid, inv_year, state_cd=23, elev=1000, slope=5, site_idx=42):
    return pd.DataFrame([{
        'stand_id':sid,'variant':'NE','inv_year':int(inv_year),
        'latitude':46.5,'longitude':-68.7,'region':9,'forest':0,'district':0,
        'basal_area_factor':0.0,'inv_plot_size':1.0,'brk_dbh':999.0,'num_plots':1,
        'age':60,'aspect':0,'slope':slope,'elevft':elev,
        'site_species':12,'site_index':site_idx,
        'state':state_cd,'county':0,'forest_type':121,'sam_wt':1.0,
    }])

def build_treeinit(t1_trees, sid):
    rows=[]
    for _,r in t1_trees.iterrows():
        try: spcd=int(r['SPCD']); dbh=float(r['DIA'])
        except: continue
        if not (math.isfinite(dbh) and dbh>=5.0): continue
        tpa=float(r['TPA_UNADJ']) if pd.notna(r.get('TPA_UNADJ')) else 6.018
        if tpa<=0: tpa=6.018
        ht=float(r['HT']) if (pd.notna(r.get('HT')) and float(r['HT'])>0) else 0.0
        cr=float(r['CR']) if pd.notna(r.get('CR')) else 50.0
        if cr<=9: cr*=10
        cr=int(max(15,min(99,cr)))
        rows.append({'stand_id':sid,'plot_id':1,'tree_id':len(rows)+1,
                     'tree_count':round(tpa,4),'species':spcd,
                     'diameter':round(dbh,2),'ht':round(ht,1),'crratio':cr})
    return pd.DataFrame(rows)

def run_fvs(sid, stand_df, tree_df, ncyc):
    binary=os.path.join(LIB_DIR,'FVSne')
    if not os.path.exists(binary): return None, f'binary not found: {binary}'
    tmp=tempfile.mkdtemp(prefix='yr0_')
    try:
        db=os.path.join(tmp,'FVS_Data.db')
        con=sqlite3.connect(db)
        stand_df.to_sql('fvs_standinit',con,if_exists='replace',index=False)
        tree_df.to_sql('fvs_treeinit',  con,if_exists='replace',index=False)
        con.close()
        key=os.path.join(tmp,'run.key')
        with open(key,'w') as fh: fh.write(KEYFILE_TMPL.format(sid=sid,db=db,ncyc=ncyc))
        subprocess.run([binary,f'--keywordfile={key}'],cwd=tmp,
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=60)
        con=sqlite3.connect(db)
        df=pd.read_sql_query('SELECT * FROM FVS_Summary2',con)
        con.close()
        return df, None
    except Exception as e: return None, str(e)
    finally: shutil.rmtree(tmp,ignore_errors=True)

# Load scorecard
sc=pd.read_csv('/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard/scorecard_vol_all.csv')
ne_def=sc[(sc['variant']=='ne')&(sc['config']=='default')&sc['MCuFt_OBS'].notna()].copy()
ne_def['vol_bias_pct']=100*(ne_def['MCuFt_PRED']/ne_def['MCuFt_OBS']-1)
ne_sorted=ne_def.sort_values('BA_OBS').reset_index(drop=True)
n=len(ne_sorted)
picks=[ne_sorted.iloc[5],ne_sorted.iloc[n//4],ne_sorted.iloc[n//2],
       ne_sorted.iloc[3*n//4],ne_sorted.iloc[n-5]]

print("=== SELECTED PLOTS ===")
for p in picks:
    print(f"  PLOT={p['PLOT']}  ST={p['STATE']}  YR1={p['MEASYEAR1']}  YR2={p['MEASYEAR2']}  "
          f"PERIOD={p['PERIOD_YR']}  BA_OBS={p['BA_OBS']:.1f}  "
          f"MCuFt_OBS={p['MCuFt_OBS']:.1f}  MCuFt_PRED={p['MCuFt_PRED']:.1f}  "
          f"bias%={p['vol_bias_pct']:+.1f}%")

TREE_COLS=['PLT_CN','STATUSCD','SPCD','DIA','HT','CR','TPA_UNADJ','VOLCFNET','VOLBFNET']

results=[]
for p in picks:
    t1_cn=str(int(float(p['PLOT'])))
    state=int(p['STATE']); abbr=STATE_ABBR.get(state)
    myr1=int(p['MEASYEAR1']); myr2=int(p['MEASYEAR2']); ncyc=int(p['PERIOD_YR'])
    if not abbr: print(f"No abbr for state {state}"); continue

    plot_f=os.path.join(FIA_DIR,f'{abbr}_PLOT.csv')
    tree_f=os.path.join(FIA_DIR,f'{abbr}_TREE.csv')
    if not os.path.exists(tree_f): print(f"{tree_f} not found"); continue

    plot_df=pd.read_csv(plot_f,usecols=['CN','COUNTYCD','PLOT','MEASYEAR'],low_memory=False)
    plot_df['CN']=plot_df['CN'].astype(str)
    m1r=plot_df[plot_df['CN']==t1_cn]
    if len(m1r)==0: print(f"T1 CN {t1_cn} not in {abbr} PLOT"); continue
    county=m1r.iloc[0]['COUNTYCD']; plotno=m1r.iloc[0]['PLOT']
    m2r=plot_df[(plot_df['COUNTYCD']==county)&(plot_df['PLOT']==plotno)&(plot_df['MEASYEAR']==myr2)]
    if len(m2r)==0: print(f"T2 not found for {t1_cn}"); continue
    t2_cn=m2r.iloc[0]['CN']

    chunks=[]
    for ch in pd.read_csv(tree_f,usecols=lambda c: c in TREE_COLS,chunksize=300_000,low_memory=False):
        ch['PLT_CN']=ch['PLT_CN'].astype(str)
        sub=ch[ch['PLT_CN'].isin({t1_cn,t2_cn})]
        if len(sub): chunks.append(sub)
    if not chunks: print(f"No trees for {t1_cn}"); continue
    all_tr=pd.concat(chunks,ignore_index=True)

    t1t=all_tr[(all_tr['PLT_CN']==t1_cn)&(all_tr['STATUSCD']==1)].copy()
    t2t=all_tr[(all_tr['PLT_CN']==t2_cn)&(all_tr['STATUSCD']==1)].copy()

    for col in ['VOLCFNET','VOLBFNET','TPA_UNADJ','DIA']:
        t1t[col]=pd.to_numeric(t1t[col],errors='coerce').fillna(0)
        if col in t2t: t2t[col]=pd.to_numeric(t2t[col],errors='coerce').fillna(0)

    fia_t1_mcuft=(t1t['VOLCFNET']*t1t['TPA_UNADJ']).sum()
    fia_t1_bdft =(t1t['VOLBFNET']*t1t['TPA_UNADJ']).sum()
    fia_t1_ba   =(t1t['TPA_UNADJ']*0.005454154*t1t['DIA']**2).sum()
    fia_t2_mcuft=(t2t['VOLCFNET']*t2t['TPA_UNADJ']).sum()

    sid=f"NE_{t1_cn[-7:]}"
    stand_df=build_standinit(sid,myr1,state_cd=state)
    tree_df=build_treeinit(t1t,sid)

    sum0,err0=run_fvs(sid,stand_df,tree_df,0)
    if err0: print(f"  FVS yr0 error {t1_cn}: {err0}"); fvs_yr0_mcuft=fvs_yr0_bdft=fvs_yr0_ba=float('nan')
    else:
        r0=sum0.iloc[0]
        fvs_yr0_mcuft=float(r0.get('MCuFt',float('nan')))
        fvs_yr0_bdft =float(r0.get('BdFt', float('nan')))
        fvs_yr0_ba   =float(r0.get('BA',   float('nan')))

    sumN,errN=run_fvs(sid,stand_df,tree_df,ncyc)
    fvs_yrN_mcuft=float('nan')
    if not errN and sumN is not None and len(sumN)>0:
        fvs_yrN_mcuft=float(sumN.iloc[-1].get('MCuFt',float('nan')))

    yr0_ratio=fvs_yr0_mcuft/fia_t1_mcuft if fia_t1_mcuft>0 else float('nan')
    yrN_ratio=fvs_yrN_mcuft/fia_t2_mcuft if fia_t2_mcuft>0 else float('nan')

    results.append({'t1_cn':t1_cn,'state':abbr,'ncyc':ncyc,
                    'fia_t1_ba':fia_t1_ba,'fvs_yr0_ba':fvs_yr0_ba,
                    'fia_t1_mcuft':fia_t1_mcuft,'fvs_yr0_mcuft':fvs_yr0_mcuft,
                    'fia_t1_bdft':fia_t1_bdft,'fvs_yr0_bdft':fvs_yr0_bdft,
                    'yr0_ratio':yr0_ratio,
                    'fia_t2_mcuft':fia_t2_mcuft,'fvs_yrN_mcuft':fvs_yrN_mcuft,
                    'yrN_ratio':yrN_ratio,
                    'sc_pred':float(p['MCuFt_PRED']),'sc_obs':float(p['MCuFt_OBS'])})

    print(f"\nPlot {t1_cn} ({abbr}, {ncyc}yr, {len(tree_df)} trees->FVS of {len(t1t)} FIA live trees):")
    print(f"  FIA T1:  BA={fia_t1_ba:.1f}  MCuFt={fia_t1_mcuft:.1f}  BdFt={fia_t1_bdft:.1f}")
    print(f"  FVS yr0: BA={fvs_yr0_ba:.1f}  MCuFt={fvs_yr0_mcuft:.1f}  BdFt={fvs_yr0_bdft:.1f}  [ratio={yr0_ratio:.3f}]")
    print(f"  FIA T2 obs:   MCuFt={fia_t2_mcuft:.1f}")
    print(f"  FVS yr{ncyc}:  MCuFt={fvs_yrN_mcuft:.1f}  [scorecard: pred={p['MCuFt_PRED']:.1f} obs={p['MCuFt_OBS']:.1f}]")
    print(f"  FVS/FIA at T2: {yrN_ratio:.3f}")
    sys.stdout.flush()

print("\n"+"="*75)
print("  YEAR-0 DIAGNOSTIC SUMMARY")
print("="*75)
print(f"  {'Plot':>22} {'St':>3} {'FIA_T1_MCuFt':>13} {'FVS_yr0':>10} "
      f"{'yr0_ratio':>10} {'FIA_T2':>10} {'FVS_yrN':>10} {'yrN_ratio':>10}")
print("-"*75)
for r in results:
    print(f"  {r['t1_cn']:>22} {r['state']:>3} {r['fia_t1_mcuft']:>13.1f} "
          f"{r['fvs_yr0_mcuft']:>10.1f} {r['yr0_ratio']:>10.3f} "
          f"{r['fia_t2_mcuft']:>10.1f} {r['fvs_yrN_mcuft']:>10.1f} {r['yrN_ratio']:>10.3f}")
if results:
    mr0=np.nanmean([r['yr0_ratio'] for r in results])
    mrN=np.nanmean([r['yrN_ratio'] for r in results])
    print(f"\n  Mean yr0 ratio (FVS/FIA): {mr0:.3f}  ({100*(mr0-1):+.1f}%)")
    print(f"  Mean yrN ratio (FVS/FIA): {mrN:.3f}  ({100*(mrN-1):+.1f}%)")
    if mr0>1.10:
        print("\n  DIAGNOSIS: Volume gap is present at YEAR 0 (initialization).")
        print("  The FVS volume equations compute more volume than FIA VOLCFNET for the same trees.")
        print("  This is NOT a growth model issue -- it is a VOLUME EQUATION mismatch.")
    elif mr0<0.90:
        print("\n  DIAGNOSIS: FVS underpredicts volume at year 0.")
    else:
        print("\n  DIAGNOSIS: Year-0 is close (<10% bias). Gap builds during projection.")
