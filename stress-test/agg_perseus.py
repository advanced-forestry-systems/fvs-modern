import glob, pandas as pd
D="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_perseus_wo1"
fs=sorted(glob.glob(D+"/perseus_100yr_agb_batch*.csv"))
print("files:",len(fs))
df=pd.concat((pd.read_csv(f) for f in fs), ignore_index=True)
print("rows:",len(df),"cols:",list(df.columns))
df.to_csv(D+"/perseus_100yr_agb_all.csv", index=False)
s=(df.groupby(["YEAR","VARIANT","CONFIG"])
     .agg(n_plots=("PLOT","nunique"), mean_agb_tons_ac=("AGB_TONS_AC","mean"))
     .reset_index().sort_values(["VARIANT","CONFIG","YEAR"]))
s.to_csv(D+"/perseus_100yr_agb_summary.csv", index=False)
print("summary rows:",len(s))
print(s.to_string(index=False))
