# -*- coding: utf-8 -*-
"""269_supp8_TS_MH2.py — TS as a weak source and MH2 dual programme."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
ecm=lr[lr.axis=="ECM-integrin"]
src_sum=ecm.groupby("source")["prob"].sum().reindex(["TS","MH2","MyoEpi","M"])
fig=plt.figure(figsize=(7.2,3.0)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.2],wspace=0.4,left=0.08,right=0.97,top=0.88,bottom=0.2)
ax=fig.add_subplot(gs[0]); ax.bar(src_sum.index,src_sum.values,color=["#4C78A8","#F58518","#54A24B","#B279A2"],width=0.62)
for i,v in enumerate(src_sum.values): ax.text(i,v+0.02,f"{v:.2f}",ha="center",fontsize=6)
ax.set_ylabel("Summed ECM-integrin probability"); ax.set_title("Source strength"); ax.set_ylim(0,src_sum.max()*1.25)
for sp in ["top","right"]: ax.spines[sp].set_visible(False)
letter(ax,"a")
# MH2 module scores
mods=pd.read_csv(os.path.join(OUT,"Fig5_coexpression_modules.csv"))
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
mh2=a.obs.cell_type=="MH2"; B=a[mh2].copy(); B.X=B.layers["lognorm"].copy()
X=np.asarray(B.X.todense(),dtype=np.float32); Z=(X-X.mean(0))/(X.std(0)+1e-9)
rows=[]
for m,lab in [(2,"M2"),(5,"M5"),(7,"M7")]:
    gids=[str(g).replace("ncbi_","") for g in mods[mods.module==m].gid]
    idx=[list(B.var_names).index(g2v[g]) for g in gids if g in g2v]
    sc=Z[:,idx].mean(1) if idx else np.zeros(B.n_obs)
    for grp in ["blank","treated"]:
        mm=(B.obs.group==grp).values
        rows.append({"module":lab,"group":grp,"mean":float(sc[mm].mean()),"n":int(mm.sum())})
d=pd.DataFrame(rows); d.to_csv(os.path.join(DEST,"Supp8_MH2_modules.csv"),index=False,encoding="utf-8-sig")
ax=fig.add_subplot(gs[1]); x=np.arange(3); w=0.34
for k,(grp,col,lab) in enumerate([("blank","#4C78A8","Control"),("treated","#E45756","Hemospasia")]):
    v=d[d.group==grp].set_index("module").reindex(["M2","M5","M7"])["mean"].values
    ax.bar(x+(k-0.5)*w,v,width=w,color=col,label=lab)
ax.set_xticks(x); ax.set_xticklabels(["M2 ECM/lipid","M5 ECM/contractile","M7 PH progenitor"],rotation=20)
ax.set_ylabel("Mean module score (MH2)"); ax.axhline(0,color="#999",lw=0.5)
ax.legend(frameon=False,fontsize=6); ax.set_title("MH2 programme balance")
for sp in ["top","right"]: ax.spines[sp].set_visible(False)
letter(ax,"b")
fig.savefig(os.path.join(DEST,"png","Supp8_TS_MH2.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(DEST,"svg","Supp8_TS_MH2.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"tiff","Supp8_TS_MH2.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
print("done")
