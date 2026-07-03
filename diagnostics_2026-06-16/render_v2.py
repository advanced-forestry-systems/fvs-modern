import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Figure 1: COND-undisturbed structural signature (BA / TPH / QMD per variant) ----
d=pd.read_csv("mv2b.csv"); d=d[d.undist_n>=15].copy()
region={"ne":"Northeast","acd":"Acadian","sn":"Southern","ls":"Lake States","cs":"Central States",
        "ie":"Inland Empire","kt":"Kootenai/Kan.","ci":"Central Idaho","tt":"Tetons","ut":"Utah",
        "ca":"Inland CA","ws":"W. Sierra","nc":"Klamath","so":"S.OR/NE CA","ec":"East Cascades",
        "wc":"West Cascades","oc":"Oregon Coast","op":"Olympic","pn":"PNW Coast"}
d["label"]=[f"{region.get(v,v)} ({v.upper()})" for v in d.variant]
d=d.sort_values("undist_QMD%").reset_index(drop=True); y=range(len(d))
fig,ax=plt.subplots(figsize=(10.2,8.4))
ax.axvspan(-5,5,color="#888888",alpha=0.08,zorder=0); ax.axvline(0,color="#888888",lw=0.7,zorder=1)
ax.scatter(d["undist_TPH%"],y,color="#d95f02",s=58,marker="v",zorder=4,label="TPH (stem density)")
ax.scatter(d["undist_BA%"],y,color="#1f78b4",s=58,marker="o",zorder=4,label="Basal area")
ax.scatter(d["undist_QMD%"],y,color="#1b7837",s=58,marker="^",zorder=4,label="QMD (mean tree size)")
for i,r in d.iterrows(): ax.plot([min(r["undist_TPH%"],r["undist_QMD%"]),max(r["undist_TPH%"],r["undist_QMD%"])],[i,i],color="#e5e5e5",lw=1.2,zorder=2)
ax.set_yticks(list(y)); ax.set_yticklabels(d["label"],fontsize=9.2)
ax.set_xticks(range(-50,40,10)); ax.set_xticklabels([f"+{x}%" if x>0 else f"{x}%" for x in range(-50,40,10)])
ax.set_xlabel("FVS bias vs observed, COND-undisturbed FIA plots only",fontsize=11)
ax.set_title("On truly undisturbed plots, FVS carries too few, too-large trees",fontsize=13,fontweight="bold",loc="left",pad=30)
ax.text(0,1.012,"Across 19 variants the signature is near-universal: QMD over-predicted (median +10%, 18/19), TPH under-predicted\n(median -14%, 15/19), basal area modestly over (+7%). Points to under-recruitment + slightly fast diameter growth.",
        transform=ax.transAxes,fontsize=9.0,color="#333",va="bottom")
ax.legend(loc="lower left",frameon=False,fontsize=9.6)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(axis="x",color="#f0f0f0",lw=0.6); ax.set_axisbelow(True)
plt.figtext(0.012,0.004,"Default FVS, one FIA remeasurement cycle. Undisturbed = FIA COND no treatment (TRTCD!=10) and no disturbance (DSTRBCD=0). Source: FIA, 2026-06-17.",fontsize=7.2,color="grey")
plt.tight_layout(rect=[0,0.018,1,0.99]); plt.savefig("cond_undisturbed_signature.png",dpi=200,facecolor="white"); print("fig1 ok")

# ---- Figure 2: BAIMULT calibration sweep on undisturbed plots ----
c=pd.read_csv("cu.csv"); lvl={"default":1.0,"BAIMULT0.90":0.90,"BAIMULT0.80":0.80,"BAIMULT0.70":0.70}
c["mult"]=c.arm.map(lvl)
vars_=["ne","sn","kt","nc"]; titles={"ne":"Northeast (NE)","sn":"Southern (SN)","kt":"Kootenai/Kan. (KT)","nc":"Klamath (NC)"}
fig2,axs=plt.subplots(1,4,figsize=(13.6,3.8),sharex=True)
for ax,v in zip(axs,vars_):
    s=c[c.variant==v].sort_values("mult",ascending=False)
    ax.axhline(0,color="#888",lw=0.7)
    ax.plot(s["mult"],s["BA%"],"-o",color="#1f78b4",label="BA",ms=5)
    ax.plot(s["mult"],s["QMD%"],"-^",color="#1b7837",label="QMD",ms=5)
    ax.plot(s["mult"],s["TPH%"],"-v",color="#d95f02",label="TPH",ms=5)
    ax.set_title(titles[v],fontsize=10.5,fontweight="bold"); ax.set_xticks([0.7,0.8,0.9,1.0]); ax.invert_xaxis()
    ax.set_xlabel("BAIMULT (diameter-growth multiplier)",fontsize=8.6)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.grid(axis="y",color="#f0f0f0",lw=0.6); ax.set_axisbelow(True)
    ax.set_yticklabels([f"+{int(t)}%" if t>0 else f"{int(t)}%" for t in ax.get_yticks()],fontsize=8)
axs[0].set_ylabel("bias vs observed",fontsize=9.5); axs[0].legend(frameon=False,fontsize=8.8,loc="center left")
fig2.suptitle("Diameter-growth calibration (BAIMULT) pulls BA and QMD toward zero; TPH (density) is untouched",fontsize=12.5,fontweight="bold",x=0.012,ha="left",y=0.99)
plt.figtext(0.012,0.005,"COND-undisturbed plots, one cycle. The growth lever fixes size (BA/QMD); the persistent TPH gap is a recruitment/ingrowth problem (a separate lever). Source: FIA, 2026-06-17.",fontsize=7.4,color="grey")
plt.tight_layout(rect=[0,0.03,1,0.94]); plt.savefig("baimult_calibration_undisturbed.png",dpi=200,facecolor="white"); print("fig2 ok")
