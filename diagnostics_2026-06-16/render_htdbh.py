import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_csv("htdbh.csv")
order=["1-3 in","3-5 in","5-9 in","9-13 in","13-19 in","19-40 in"]
d["xi"]=d.dbh_class.map({c:i for i,c in enumerate(order)})
region={"ne":"Northeast (NE)","sn":"Southern (SN)","kt":"Kootenai (KT)","ie":"Inland Empire (IE)","pn":"PNW Coast (PN)"}
col={"ne":"#1b7837","sn":"#d95f02","kt":"#7570b3","ie":"#b2182b","pn":"#1f78b4"}
fig,ax=plt.subplots(figsize=(10.4,6.2))
ax.axhspan(-5,5,color="#888",alpha=0.08,zorder=0); ax.axhline(0,color="#666",lw=0.8,zorder=1)
for v in ["ne","sn","kt","ie","pn"]:
    s=d[d.variant==v].sort_values("xi")
    if len(s)==0: continue
    ax.plot(s.xi,s["ht_bias%"],"-o",color=col[v],lw=2,ms=6,label=region[v],zorder=4)
# pooled (from full run)
pooled={0:-6.8,1:1.8,2:3.7,3:3.8,4:6.9,5:8.8}
ax.plot(list(pooled.keys()),list(pooled.values()),"--",color="#222",lw=2.4,label="All-variant pooled",zorder=5)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
ax.set_xlabel("DBH class",fontsize=11); ax.set_ylabel("FVS HT-DBH height bias vs FIA measured (%)",fontsize=11)
ax.set_title("FVS HT-DBH curve bias is variant- and size-specific (hidden by unbiased top height)",fontsize=12.5,fontweight="bold",loc="left",pad=12)
ax.legend(frameon=False,fontsize=9.5,ncol=2,loc="upper left")
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(axis="y",color="#f0f0f0",lw=0.6); ax.set_axisbelow(True)
ax.set_yticklabels([f"+{int(t)}%" if t>0 else f"{int(t)}%" for t in ax.get_yticks()])
plt.figtext(0.012,0.005,"FVS heights imputed via HT-DBH (input heights blanked) vs FIA field-measured ACTUALHT, by DBH class. IE/KT over-predict; SN under-predicts large trees; PN curve shape wrong; NE mild. 82 per-species recal ratios derived. Source: FIA, 2026-06-17.",fontsize=7.2,color="grey")
plt.tight_layout(rect=[0,0.02,1,1]); plt.savefig("htdbh_curve_bias.png",dpi=200,facecolor="white"); print("ok")
