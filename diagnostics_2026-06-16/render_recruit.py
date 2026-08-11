import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
dec=pd.read_csv("dec.csv"); bench=pd.read_csv("mv2b.csv")[["variant","undist_TPH%"]]
d=dec.merge(bench,on="variant",how="left")
region={"ne":"NE","acd":"ACD","sn":"SN","ls":"LS","kt":"KT","ci":"CI","nc":"NC","ec":"EC","pn":"PN"}
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12.6,5.4))

# Left: observed ingrowth rate vs FVS undisturbed TPH bias
x=d["ingrowth_rate_pct_per_decade"]; y=d["undist_TPH%"]
ax1.scatter(x,y,s=90,color="#d95f02",zorder=4,edgecolor="white",linewidth=0.8)
for _,r in d.iterrows():
    ax1.annotate(region.get(r.variant,r.variant),(r["ingrowth_rate_pct_per_decade"],r["undist_TPH%"]),
                 xytext=(4,4),textcoords="offset points",fontsize=9,color="#333")
m=(~x.isna())&(~y.isna()); b1,b0=np.polyfit(x[m],y[m],1); xs=np.linspace(x.min(),x.max(),50)
ax1.plot(xs,b0+b1*xs,"--",color="#762a83",lw=1.6,zorder=3)
r2=np.corrcoef(x[m],y[m])[0,1]**2
ax1.axhline(0,color="#aaa",lw=0.6)
ax1.set_xlabel("Observed ingrowth rate (% of initial TPH per decade)",fontsize=10.5)
ax1.set_ylabel("FVS undisturbed TPH bias (%)",fontsize=10.5)
ax1.set_title("More real ingrowth -> worse FVS density",fontsize=12,fontweight="bold",loc="left")
ax1.text(0.97,0.06,"R$^2$=%.2f"%r2,transform=ax1.transAxes,ha="right",fontsize=11,color="#762a83")
for s in ["top","right"]: ax1.spines[s].set_visible(False)
ax1.grid(color="#f0f0f0",lw=0.6); ax1.set_axisbelow(True)

# Right: net TPH change observed (ingrowth - mortality) vs FVS (mortality only)
vv=d.sort_values("ingrowth_rate_pct_per_decade"); yy=range(len(vv))
ax2.axvline(0,color="#888",lw=0.7)
ax2.barh([i+0.2 for i in yy],vv["obs_net_TPH"],height=0.38,color="#1b7837",label="Observed net (ingrowth - mortality)")
ax2.barh([i-0.2 for i in yy],vv["fvs_net_TPH"],height=0.38,color="#b2182b",label="FVS net (mortality only)")
ax2.set_yticks(list(yy)); ax2.set_yticklabels([region.get(v,v) for v in vv.variant],fontsize=9.5)
ax2.set_xlabel("Net change in stems/ha over the interval",fontsize=10.5)
ax2.set_title("Real stands self-replace; FVS only loses stems",fontsize=12,fontweight="bold",loc="left")
ax2.legend(frameon=False,fontsize=9,loc="lower right")
for s in ["top","right"]: ax2.spines[s].set_visible(False)
ax2.grid(axis="x",color="#f0f0f0",lw=0.6); ax2.set_axisbelow(True)

fig.suptitle("The undisturbed TPH gap is missing recruitment: FVS has no background ingrowth in undisturbed stands",
             fontsize=13,fontweight="bold",x=0.012,ha="left",y=0.99)
plt.figtext(0.012,0.005,"COND-undisturbed FIA plots, one remeasurement cycle, default FVS. FVS establishment model is disturbance-triggered (partial model), so undisturbed runs add no ingrowth. Source: FIA, 2026-06-17.",fontsize=7.4,color="grey")
plt.tight_layout(rect=[0,0.03,1,0.95]); plt.savefig("recruitment_gap.png",dpi=200,facecolor="white"); print("ok r2=%.2f"%r2)
