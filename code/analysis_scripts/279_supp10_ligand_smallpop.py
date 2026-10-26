# -*- coding: utf-8 -*-
"""279_supp10_ligand_smallpop.py — VWF/hemocytin and THBS4 in small niche populations."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad")); T=a.obsm["X_tsne"]
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(gid): return np.asarray(a[:,g2v[gid]].layers["lognorm"].todense()).ravel()
vwf=expr("113824171"); thbs4=expr("113823118")
fig,axes=plt.subplots(2,2,figsize=(6.8,5.2))
for k,(v,title) in enumerate([(vwf,"VWF/hemocytin"),(thbs4,"THBS4")]):
    ax=axes[k,0]
    vmax=np.quantile(v,0.99) or 1
    sc=ax.scatter(T[:,0],T[:,1],c=v,cmap="Reds",s=1.0,linewidths=0,rasterized=True,vmin=0,vmax=vmax)
    for ct in ["MyoEpi","M","MH2","TS"]:
        m=(a.obs.cell_type==ct).values
        ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=5.8,ha="center",bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title+" (t-SNE)")
    for s in ax.spines.values(): s.set_visible(False)
    plt.colorbar(sc,ax=ax,fraction=0.045,pad=0.02).set_label("log1p CP10K",fontsize=6)
    letter(ax,"ac"[k])
for k,(v,title) in enumerate([(vwf,"VWF/hemocytin"),(thbs4,"THBS4")]):
    ax=axes[k,1]
    cts=["MyoEpi","M","MH2","TS"]; data=[v[(a.obs.cell_type==ct).values] for ct in cts]
    parts=ax.violinplot(data,showextrema=False,widths=0.8)
    for i,b in enumerate(parts["bodies"]): b.set_facecolor("#C44E52" if i==0 else "#8DA0CB"); b.set_alpha(0.75); b.set_edgecolor("#444"); b.set_linewidth(0.4)
    for i,d in enumerate(data): ax.plot([i,i],[np.median(d),np.median(d)],color="#222",lw=1.0)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts,rotation=30); ax.set_ylabel("log1p CP10K"); ax.set_title(title+" by source")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    letter(ax,"bd"[k])
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,ext,"Supp10_ligand_smallpop."+ext),bbox_inches="tight",**kw)
print("saved")
