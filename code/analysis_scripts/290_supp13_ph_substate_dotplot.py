# -*- coding: utf-8 -*-
"""290_supp13_ph_substate_dotplot.py — PH substate marker dot plot + naming audit."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
a=ad.read_h5ad(os.path.join(OUT,"step3","step3_lyo_scores.h5ad"),backed="r")
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv")); ph.columns=["cell","group","sub","u1","u2"]
sub=ph.set_index("cell")["sub"]; keep=np.isin(a.obs_names,sub.index)
A=a[keep]; s=sub.reindex(A.obs_names).values
genes=[("s0","PEN-3","ncbi_113808997"),("s0","LOC113819543","ncbi_113819543"),("s0","LOC113819837","ncbi_113819837"),("s0","LOC113830563","ncbi_113830563"),("s0","LOC113802679","ncbi_113802679"),
       ("s1","LOC113808957","ncbi_113808957"),("s1","LOC113808628","ncbi_113808628"),("s1","LOC113800422","ncbi_113800422"),("s1","Igfbp7","ncbi_113811702"),("s1","ccnd2","ncbi_113814652"),
       ("s2","LOC113810024","ncbi_113810024"),("s2","THBS4","ncbi_113823118"),("s2","TGM1","ncbi_113823945"),
       ("s3","gene-LOC113818748","gene-LOC113818748"),("s3","gene-LOC113823090","gene-LOC113823090"),("s3","gene-LOC113823089","gene-LOC113823089"),
       ("context","PPAF2","ncbi_113812814"),("context","Grik2","ncbi_113823802")]
rows=[]
for grp,name,key in genes:
    hit=None
    if key in A.var_names: hit=key
    elif key.replace("gene-","") in A.var_names: hit=key.replace("gene-","")
    if hit is None: print("missing",name,key); continue
    x=np.asarray(A[:,hit].layers["lognorm"].todense()).ravel()
    for s_ in [0,1,2,3]:
        m=s==s_; rows.append({"group":grp,"gene":name,"sub":s_,"mean":x[m].mean(),"pct":(x[m]>0).mean()*100})
d=pd.DataFrame(rows); d.to_csv(os.path.join(DEST,"Supp13_PH_substate_markers.csv"),index=False,encoding="utf-8-sig")
order=d.drop_duplicates("gene").gene.tolist()
fig,ax=plt.subplots(figsize=(4.8,5.2))
X=[];Y=[];SZ=[];CC=[]
for i,g in enumerate(order):
    for j,s_ in enumerate([0,1,2,3]):
        r=d[(d.gene==g)&(d["sub"]==s_)]
        if len(r)==0: continue
        r=r.iloc[0]; X.append(j); Y.append(i); SZ.append(max(float(r["pct"])*1.1,2.0)); CC.append(float(r["mean"]))
ax.scatter(X,Y,s=SZ,c=CC,cmap="Reds",vmin=0,vmax=float(d["mean"].max()),edgecolors="#666",linewidths=0.25)
ax.set_xticks([0,1,2,3]); ax.set_xticklabels(["s0\nECM/ribosome","s1\nprogenitor","s2\nimmune-primed","s3\nTGM1-low"],fontsize=6.2)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=6)
ax.set_xlim(-0.6,3.6); ax.set_ylim(len(order)-0.5,-0.5); ax.set_title("PH substate markers")
for v in [0,25,50,75]: ax.scatter([],[],s=max(v*1.1,2),c="#BBB",edgecolors="#666",linewidths=0.25,label=f"{v}%")
ax.legend(title="% expressing",frameon=False,fontsize=5.2,title_fontsize=5.2,loc="center left",bbox_to_anchor=(1.01,0.5))
for sp in ["top","right"]: ax.spines[sp].set_visible(False)
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,ext,"Supp13_PH_substate_markers."+ext),bbox_inches="tight",**kw)
print("saved")
