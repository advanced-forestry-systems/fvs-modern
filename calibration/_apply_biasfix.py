f="/users/PUOM0008/crsfaaron/fvs-modern/calibration/run_5model_scorecard_biasfix.py"
s=open(f).read()
old='''def obs_ba_ft2ac(tree_df: pd.DataFrame, plt_cn: int) -> float:
    """Observed BA (ft2/ac) from FIA TREE records at the remeasurement year."""
    t = tree_df.loc[tree_df["PLT_CN"] == plt_cn].copy()
    t = t.loc[t["DIA"].notna() & t["TPA_UNADJ"].notna()]
    if len(t) == 0:
        return float("nan")
    return float((t["TPA_UNADJ"] * 0.005454154 * t["DIA"] ** 2).sum())'''
new='''def obs_ba_ft2ac(tree_df: pd.DataFrame, plt_cn: int, min_dia: float = 5.0) -> float:
    """Observed BA (ft2/ac) from FIA TREE records at the remeasurement year.

    BIASFIX: apply the same DIA >= min_dia (default 5.0 in) threshold used by
    build_treeinit() for the FVS input tree list. FVS is only given >=5" trees
    and its Summary BA reflects only those; the raw FIA TREE table also carries
    1-5" microplot saplings at ~75 TPA each that add basal area but carry NULL
    merch volume. Counting saplings in OBS BA (but not in PRED BA, nor in OBS
    volume) manufactured a spurious ~20-30% BA under-prediction while leaving
    volume comparisons clean -- the contradictory signature Aaron flagged.
    """
    t = tree_df.loc[tree_df["PLT_CN"] == plt_cn].copy()
    t = t.loc[t["DIA"].notna() & t["TPA_UNADJ"].notna() & (t["DIA"] >= min_dia)]
    if len(t) == 0:
        return float("nan")
    return float((t["TPA_UNADJ"] * 0.005454154 * t["DIA"] ** 2).sum())'''
assert old in s, "FUNC not found"
s=s.replace(old,new)
old2='''        ba_obs              = obs_ba_ft2ac(trees, t2_cn)
        mcuft_obs, bdft_obs = obs_vol_mcuft_bdft(trees, t2_cn)'''
new2='''        ba_obs              = obs_ba_ft2ac(trees, t2_cn)              # BIASFIX >=5"
        ba_obs_allsize      = obs_ba_ft2ac(trees, t2_cn, min_dia=0.0)  # diag: old all-size
        mcuft_obs, bdft_obs = obs_vol_mcuft_bdft(trees, t2_cn)'''
assert old2 in s, "CALL not found"
s=s.replace(old2,new2)
old3='''                "BA_OBS":     ba_obs,
                "MCuFt_OBS":  mcuft_obs,'''
new3='''                "BA_OBS":         ba_obs,
                "BA_OBS_ALLSIZE": ba_obs_allsize,
                "MCuFt_OBS":  mcuft_obs,'''
assert old3 in s, "ROW not found"
s=s.replace(old3,new3)
open(f,"w").write(s)
print("PATCHED bytes",len(s))
