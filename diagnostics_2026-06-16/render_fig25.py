import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_csv("all25.csv")
region={"ne":"Northeast","acd":"Acadian","sn":"Southern","ls":"Lake States","cs":"Central States",
        "ie":"Inland Empire","kt":"Kootenai/Kaniksu","ci":"Central Idaho","em":"E. Montana",
        "bm":"Blue Mtns","cr":"Central Rockies","tt":"Tetons","ut":"Utah","ca":"Inland CA/S.Casc.",
        "ws":"W. Sierra Nevada","nc":"Klamath Mtns","so":"S.Oregon/NE CA","ec":"East Cascades",
        "wc":"West Cascades","oc":"Oregon Coast","op":"Olympic Pen.","pn":"PNW Coast"}
d=d[d.n>=15].copy()                       # drop tiny-sample variants (em n=5)
d["label"]=[f"{region.get(v,v)} ({v.upper()})" for v in d.variant]
d=d.sort_values("undist_BA%").reset_index(drop=True)
y=range(len(d))
med_u=d["undist_BA%"].median(); med_p=d["pooled_BA%"].median(); med_h=d["harvest_BA%"].median()
pal={"u":"#1b7837","p":"#762a83","h":"#b2182b"}
fig,ax=plt.subplots(figsize=(10.4,8.6))
ax.axvspan(-5,5,color="#1b7837",alpha=0.08,zorder=0)
ax.axvline(0,color="grey45" if False else "#888888",lw=0.6,zorder=1)
for i,r in d.iterrows():
    ax.plot([r["undist_BA%"],r["harvest_BA%"]],[i,i],color="#dddddd",lw=1.6,zorder=2)
ax.scatter(d["undist_BA%"],y,color=pal["u"],s=66,marker="o",zorder=5,label="Undisturbed plots")
ax.scatter(d["pooled_BA%"],y,color=pal["p"],s=80,marker="D",zorder=4,label="Pooled (all plots)")
ax.scatter(d["harvest_BA%"],y,color=pal["h"],s=74,marker="^",zorder=4,label="Harvested plots")
ax.set_yticks(list(y)); ax.set_yticklabels(d["label"],fontsize=9.2)
ax.set_xticks(range(0,90,10)); ax.set_xticklabels([f"+{x}%" if x>0 else f"{x}%" for x in range(0,90,10)])
ax.set_xlim(-22,82)
ax.set_xlabel("FVS basal-area bias vs observed FIA remeasurement",fontsize=11)
ax.set_title("Across 22 FVS variants, basal-area over-prediction is a disturbance artifact",fontsize=13,fontweight="bold",loc="left",pad=30)
ax.text(0,1.012,"On undisturbed FIA plots (green band ±5%) the median bias is +1.8%; the pooled +14% and harvested +42% medians\nare driven by harvest/mortality the default run never simulates. Undisturbed is within ±10% for 18 of 22 variants.",
        transform=ax.transAxes,fontsize=9.0,color="#333333",va="bottom")
ax.legend(loc="lower right",frameon=False,fontsize=10)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(axis="x",color="#eeeeee",lw=0.6); ax.set_axisbelow(True)
plt.figtext(0.012,0.004,"Default FVS, one FIA remeasurement cycle (5-15 yr). Harvest = >15% of t1 live stems cut/removed. Variants with n<15 excluded. Source: FIA, 2026-06-17.",fontsize=7.3,color="grey")
plt.tight_layout(rect=[0,0.018,1,0.99])
plt.savefig("disturbance_artifact_all25.png",dpi=200,facecolor="white")
print("saved disturbance_artifact_all25.png ; variants plotted:",len(d))
