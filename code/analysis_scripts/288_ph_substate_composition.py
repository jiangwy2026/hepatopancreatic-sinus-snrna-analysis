# -*- coding: utf-8 -*-
"""288_ph_substate_composition.py — PH substate composition stacked bar."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig2","preview"); os.makedirs(DEST,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":8,
 "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7.5,"ytick.labelsize":7.5,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})

d=pd.read_csv(os.path.join(OUT,"Fig5_dynamic_subpop_subcluster_composition.csv"))
d=d[d.population=="PH"].sort_values("subcluster")
labels={0:"PH_s0\nECM/ribosome-high",1:"PH_s1\nprogenitor-like",2:"PH_s2\nimmune-primed",3:"PH_s3\nTGM1+"}
cols={0:"#4C78A8",1:"#F58518",2:"#E45756",3:"#54A24B"}
fig,ax=plt.subplots(figsize=(3.8,3.4))
groups=["Control","Hemospasia"]
bottom=np.zeros(2)
for _,r in d.iterrows():
    vals=np.array([r.pct_control,r.pct_hemospasia]); ns=np.array([r.n_control,r.n_hemospasia])
    xpos=[0,0.75]
    ax.bar(xpos,vals,bottom=bottom,color=cols[int(r.subcluster)],width=0.62,edgecolor="white",linewidth=0.8,label=labels[int(r.subcluster)].split("\n")[0])
    for i,(v,n,b) in enumerate(zip(vals,ns,bottom)):
        if v>6: ax.text(xpos[i],b+v/2,f"s{int(r.subcluster)}  {v:.1f}%\n(n={n})",ha="center",va="center",fontsize=6.5,color="white" if int(r.subcluster)!=3 else "#222")
    bottom+=vals
ax.set_xticks(xpos); ax.set_xticklabels(groups); ax.set_xlim(-0.5,1.25)
ax.set_ylim(0,100); ax.set_ylabel("PH cells (%)"); ax.set_title("PH substate composition")
ax.text(xpos[0],-8,f"n = {int(d.n_control.sum())}",ha="center",va="top",fontsize=7)
ax.text(xpos[1],-8,f"n = {int(d.n_hemospasia.sum())}  (2.30×)",ha="center",va="top",fontsize=7)
ax.legend(handles=[Patch(color=cols[k],label=labels[k].replace("\n"," ")) for k in sorted(labels)],
          frameon=False,fontsize=5.8,loc="center left",bbox_to_anchor=(1.01,0.5))
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,"Fig2_PH_substate_composition."+ext),bbox_inches="tight",**kw)
print("saved")
