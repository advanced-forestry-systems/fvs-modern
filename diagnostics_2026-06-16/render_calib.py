import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_csv("calib_full.csv")
region={"ne":"Northeast","acd":"Acadian","sn":"Southern","ls":"Lake States","cs":"Central States",
        "ie":"Inland Empire","kt":"Kootenai","ci":"Central Idaho","cr":"Central Rockies","ut":"Utah",
        "ca":"Inland CA","nc":"Klamath","ec":"East Cascades","wc":"West Cascades","pn":"PNW Coast"}
d["label"]=[region.get(v,v) for v in d.variant]
d=d.sort_values("def_QMD%").reset_index(drop=True)
metrics=[("BA","def_BA%","cal_BA%"),("TPH","def_TPH%","cal_TPH%"),("QMD","def_QMD%","cal_QMD%")]
fig,axs=plt.subplots(1,3,figsize=(14.5,7.2),sharey=True)
y=range(len(d))
for ax,(name,dc,cc) in zip(axs,metrics):
    ax.axvspan(-5,5,color="#1b7837",alpha=0.07,zorder=0); ax.axvline(0,color="#888",lw=0.6,zorder=1)
    for i,r in d.iterrows():
        ax.annotate("",xy=(r[cc],i),xytext=(r[dc],i),arrowprops=dict(arrowstyle="->",color="#bbb",lw=1.4),zorder=2)
    ax.scatter(d[dc],y,color="#b2182b",s=46,zorder=4,label="default")
    ax.scatter(d[cc],y,color="#1b7837",s=46,zorder=4,label="calibrated")
    ax.set_title(name+" bias",fontsize=12,fontweight="bold")
    ax.set_xlabel("% bias vs observed",fontsize=10)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.grid(axis="x",color="#f2f2f2",lw=0.6); ax.set_axisbelow(True)
axs[0].set_yticks(list(y)); axs[0].set_yticklabels(d["label"],fontsize=9.3); axs[0].legend(frameon=False,fontsize=9,loc="lower left")
fig.suptitle("Fully-calibrated FVS vs default across 15 variants (COND-undisturbed FIA): brms SDImax + ingrowth + BAIMULT",
             fontsize=13,fontweight="bold",x=0.012,ha="left",y=0.99)
plt.figtext(0.012,0.005,"One remeasurement cycle. Calibrated = plot-level brms max SDI + per-variant ingrowth injection + BAIMULT 0.90. Green band +/-5%. QMD/TPH improve broadly; TPH over-corrects where it was already over-predicted (sign-aware injection needed). Source: FIA, 2026-06-17.",fontsize=7.2,color="grey")
plt.tight_layout(rect=[0,0.02,1,0.96]); plt.savefig("calibrated_vs_default_allvariants.png",dpi=200,facecolor="white"); print("ok")
