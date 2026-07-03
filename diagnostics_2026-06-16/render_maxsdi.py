import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fia=pd.read_csv("maxsdi_fia.csv"); lng=pd.read_csv("maxsdi_long.csv")
region={"ne":"Northeast","acd":"Acadian","sn":"Southern","ls":"Lake States","cs":"Central States",
        "ie":"Inland Empire","kt":"Kootenai","ci":"Central Idaho","cr":"Central Rockies","ut":"Utah",
        "ca":"Inland CA","nc":"Klamath","ec":"East Cascades","wc":"West Cascades","pn":"PNW Coast"}
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13.4,6.0))

# Panel A: FIA observed max SDI per variant (p50-p99 range, p95 marker)
f=fia.sort_values("SDI_p95").reset_index(drop=True); y=range(len(f))
for i,r in f.iterrows():
    ax1.plot([r.SDI_p50,r.SDI_p99],[i,i],color="#cfe3d6",lw=5,solid_capstyle="round",zorder=2)
ax1.scatter(f.SDI_p95,y,color="#1b7837",s=55,zorder=4,label="p95 (proposed max SDI)")
ax1.scatter(f.SDI_p50,y,color="#9ecae1",s=28,zorder=3,label="p50 (median stocking)")
ax1.scatter(f.SDI_p99,y,color="#762a83",s=28,marker="|",zorder=3,label="p99")
ax1.set_yticks(list(y)); ax1.set_yticklabels([region.get(v,v) for v in f.variant],fontsize=9.2)
ax1.set_xlabel("Reineke SDI (English, per acre)",fontsize=10.5)
ax1.set_title("FIA-observed max SDI by variant",fontsize=12,fontweight="bold",loc="left")
ax1.legend(frameon=False,fontsize=8.6,loc="lower right")
for s in ["top","right"]: ax1.spines[s].set_visible(False)
ax1.grid(axis="x",color="#f0f0f0",lw=0.6); ax1.set_axisbelow(True)

# Panel B: long-term leverage - default vs revised 100-yr BA
l=lng.sort_values("BA_reduction_pct").reset_index(drop=True); yy=np.arange(len(l))
ax2.barh(yy+0.2,l.def_100yr_BA_m2ha,height=0.38,color="#b2182b",label="Default SDIMAX")
ax2.barh(yy-0.2,l.rev_100yr_BA_m2ha,height=0.38,color="#1b7837",label="Revised SDIMAX = FIA p95")
for i,r in l.iterrows():
    ax2.annotate("%.0f%%"%r.BA_reduction_pct,(r.def_100yr_BA_m2ha,i),xytext=(4,0),textcoords="offset points",va="center",fontsize=8.4,color="#b2182b")
ax2.set_yticks(list(yy)); ax2.set_yticklabels([region.get(v,v) for v in l.variant],fontsize=9.5)
ax2.set_xlabel("100-year stand basal area (m$^2$/ha)",fontsize=10.5)
ax2.set_title("Revised max SDI cuts long-term basal area 5-37%",fontsize=12,fontweight="bold",loc="left")
ax2.legend(frameon=False,fontsize=9,loc="lower right")
for s in ["top","right"]: ax2.spines[s].set_visible(False)
ax2.grid(axis="x",color="#f0f0f0",lw=0.6); ax2.set_axisbelow(True)

fig.suptitle("Max SDI is the dominant long-term control: FVS defaults exceed the FIA self-thinning limit in several variants",
             fontsize=13,fontweight="bold",x=0.012,ha="left",y=0.99)
plt.figtext(0.012,0.005,"Left: Reineke summation SDI on undisturbed FIA plots (n shown in data). Right: dense undisturbed stands projected 100 yr (10x10) default vs SDIMAX set to the FIA p95. Source: FIA, 2026-06-17.",fontsize=7.3,color="grey")
plt.tight_layout(rect=[0,0.03,1,0.95]); plt.savefig("maxsdi_longterm.png",dpi=200,facecolor="white"); print("ok")
