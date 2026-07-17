src="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/build_singlecohort_traj_v5.R"
dst="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/build_singlecohort_traj_v6.R"
s=open(src).read()
anchor='v5d  <- file.path(EQ,"out_v5_seed")'
assert s.count(anchor)==1, s.count(anchor)
s=s.replace(anchor, anchor+'\nv6d  <- file.path(EQ,"out_v6_seed")',1)
v5line='  v5_b2             = list(metrics=lf(v5d,"_conus_b2_metrics\\\\.csv$"), tl=NULL, cfg=NULL, ht_col="HT_M_MEAN"),'
assert s.count(v5line)==1, ("V5LINE", s.count(v5line))
v6lines=(v5line
  +'\n  v6_b2             = list(metrics=lf(v6d,"_conus_b2_metrics\\\\.csv$"), tl=NULL, cfg=NULL, ht_col="HT_M_MEAN"),'
  +'\n  v6_b2_dom         = list(metrics=lf(v6d,"_conus_b2_metrics\\\\.csv$"), tl=NULL, cfg=NULL, ht_col="HT_M_DOM"),')
s=s.replace(v5line, v6lines,1)
open(dst,"w").write(s)
print("WROTE",dst,"len",len(s))
print("v6d=",s.count('v6d  <- file.path'),"v6_b2 reg=",s.count('  v6_b2 '),"v6_b2_dom reg=",s.count('  v6_b2_dom '))
