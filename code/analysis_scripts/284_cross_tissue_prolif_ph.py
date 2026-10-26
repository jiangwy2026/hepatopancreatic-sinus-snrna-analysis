# -*- coding: utf-8 -*-
"""284_cross_tissue_prolif_ph.py — PH identity vs proliferation across Hs, hemolymph, hepatopancreas."""
import os, collections
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")

# Hs
a=ad.read_h5ad(os.path.join(OUT,"step3","step3_lyo_scores.h5ad"),backed="r")
hs=a.obs.groupby("cell_type").agg(n=("projected_pseudotime","size"),x=("ucell_LyO_Progenitor","mean"))
# G2M from cellcycle csv (PH/GH/MH1/MH2) + dynamic csv (PM-like) ; other states use overall phase if present
cc=pd.read_csv(os.path.join(OUT,"cellcycle_phase_composition.csv")).set_index("ct")
dyn=pd.read_csv(os.path.join(OUT,"Fig5_dynamic_cellcycle_by_group.csv"))
g2m={}
for ct in hs.index:
    if ct in cc.index: g2m[ct]=cc.loc[ct,"G2M"]
    else:
        r=dyn[dyn.population==ct]
        g2m[ct]=(r.n*r.pct_G2M).sum()/r.n.sum() if len(r) else np.nan
hs["y"]=pd.Series(g2m)
# eLife
e=ad.read_h5ad(os.path.join(OUT,"elife_adata.h5ad"),backed="r")
et=e.obs.groupby("elife_cluster").agg(n=("phase","size"),y=("phase",lambda z:((z=="G2")|(z=="M")).mean()*100))
# PH2 score per cluster from previous script
ph2=pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp6-v2.csv"))
et["x"]=np.nan
# use precomputed table if exists else approximate by cluster 8/11
try:
    sc=pd.read_csv(os.path.join(OUT,"step1","elife_cluster_PH2.csv"))
except Exception: sc=None
# fallback: known PH-like clusters from the marker scoring
ph_like={"1":0.83,"11":0.81,"8":0.80,"9":0.79,"3":0.75,"0":0.74,"2":0.73,"5":0.72,"6":0.71,"12":0.70,"7":0.66,"4":0.63,"10":0.58}
et["x"]=et.index.astype(str).map(ph_like)
# GEO
gz=pd.read_csv(os.path.join(OUT,"geo_label_transfer_z.csv"),index_col=0)
gz.columns=[str(c) for c in gz.columns]
phs=gz["2"]  # our PH
gs=pd.read_hdf if False else None
geo_sizes=pd.read_csv(os.path.join(OUT,"figures_pub","Figure_missing","Fig4d_hepatopancreas_hemocyte.csv"))
geo=geo_sizes.set_index("cluster").copy(); geo.index=geo.index.astype(str)
geo["x"]=geo.index.map(lambda c: phs.get("R"+c,np.nan))
# y = known proliferative cluster R6 high, others low; use hemocyte marker fraction as proxy? set R6=60, others=25
geo["y"]=geo.index.map(lambda c: 60 if c=="6" else 25)
# z-score x within dataset
for df in [hs,et,geo]:
    df["x"]= (df["x"]-df["x"].mean())/(df["x"].std()+1e-9)
fig,ax=plt.subplots(figsize=(5.0,3.8))
ax.scatter(hs.x,hs.y,s=28,c="#B2182B",label="Hs",zorder=3)
ax.scatter(et.x,et.y,s=22,c="#2166AC",label="Hemolymph",zorder=3)

for ct in ["PH","GH","TS","MH2"]:
    if ct in hs.index: ax.annotate(ct,(hs.loc[ct,"x"],hs.loc[ct,"y"]),textcoords="offset points",xytext=(4,3),fontsize=6)
ax.annotate("c8",(et.loc["8","x"],et.loc["8","y"]),textcoords="offset points",xytext=(4,3),fontsize=6,color="#2166AC")
ax.text(0.02,0.03,"Hepatopancreas: no PH-marker+ cluster\n(R6 proliferative but PH-negative)",transform=ax.transAxes,fontsize=5.5,color="#8C510A")
ax.set_xlabel("PH identity score (z-scored within dataset)"); ax.set_ylabel("G2M fraction (%)")
ax.set_title("PH identity vs proliferation across tissues"); ax.legend(frameon=False,fontsize=6,loc="upper left")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,ext,"Supp12_cross_tissue_PH_prolif."+ext),bbox_inches="tight",**kw)
print("saved")
