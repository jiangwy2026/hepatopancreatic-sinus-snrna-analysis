# -*- coding: utf-8 -*-
"""281_lr_sensitivity.py — permutation depth, subsampling stability, threshold sensitivity."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
labels=a.obs.cell_type.astype(str).values

# 1) permutation depth: top 6 unique pairs
top=lr.sort_values("prob",ascending=False).drop_duplicates(["pair","source","target"]).head(6)
rng=np.random.default_rng(42); rows=[]
for r in top.itertuples():
    lg=expr(str(r.ligand_genes).split(";")); rg=expr(str(r.receptor_genes).split(";"))
    src=(labels==r.source); tgt=(labels==r.target); obs=lg[src].mean()*rg[tgt].mean()
    n=5000; null=np.empty(n)
    for k in range(n):
        perm=rng.permutation(labels)
        null[k]=lg[perm==r.source].mean()*rg[perm==r.target].mean()
    p=(np.sum(null>=obs)+1)/(n+1)
    rows.append({"pair":r.pair,"source":r.source,"target":r.target,"prob":obs,"p_1000":r.p_perm,"p_5000":p})
perm=pd.DataFrame(rows); perm.to_csv(os.path.join(DEST,"Supp11_permutation_depth.csv"),index=False,encoding="utf-8-sig")

# 2) subsampling stability: 100 x n=40 per group
pairs=lr.drop_duplicates(["pair","source","target"]).copy()
vec={}
for r in pairs.itertuples():
    vec[(r.pair,r.source,r.target)]=(expr(str(r.ligand_genes).split(";")),expr(str(r.receptor_genes).split(";")))
topsrc={}
rng=np.random.default_rng(7)
for it in range(100):
    score={}
    for (p,s,t),(lg,rg) in vec.items():
        si=np.where(labels==s)[0]; ti=np.where(labels==t)[0]
        si=rng.choice(si,min(40,len(si)),replace=False); ti=rng.choice(ti,min(40,len(ti)),replace=False)
        score[(s,t)]=score.get((s,t),0)+lg[si].mean()*rg[ti].mean()
    for k,v in sorted(score.items(),key=lambda x:-x[1])[:3]: topsrc[k]=topsrc.get(k,0)+1
stab=pd.DataFrame([{"source":k[0],"target":k[1],"top3_freq":v} for k,v in topsrc.items()]).sort_values("top3_freq",ascending=False)
stab.to_csv(os.path.join(DEST,"Supp11_subsampling_stability.csv"),index=False,encoding="utf-8-sig")

# 3) threshold sensitivity
thr=[0.01,0.02,0.05,0.10,0.15]
counts=[int((lr.prob>t).sum()) for t in thr]

fig,axes=plt.subplots(1,3,figsize=(7.4,3.0),gridspec_kw={"width_ratios":[1,1.25,0.9]})
ax=axes[0]
ax.scatter(-np.log10(perm.p_1000),-np.log10(perm.p_5000),s=26,c="#C44E52",edgecolors="#333",linewidths=0.3,zorder=3)
mx=max((-np.log10(perm.p_1000)).max(),(-np.log10(perm.p_5000)).max())*1.1
ax.plot([0,mx],[0,mx],color="#999",ls="--",lw=0.6)
ax.set_xlabel("−log10 P (1,000 perms)"); ax.set_ylabel("−log10 P (5,000 perms)"); ax.set_title("Permutation depth")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"a")
ax=axes[1]
stab=stab.head(8).iloc[::-1]
ax.barh(stab.source+"→"+stab.target,stab.top3_freq,color="#4C78A8",height=0.62)
ax.set_xlabel("Top-3 frequency across 100 subsamples"); ax.set_title("Subsampling stability")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"b")
ax=axes[2]
ax.bar([str(t) for t in thr],counts,color="#54A24B",width=0.62)
for i,v in enumerate(counts): ax.text(i,v+0.5,str(v),ha="center",fontsize=6)
ax.set_xlabel("Probability threshold"); ax.set_ylabel("Passing pairs"); ax.set_title("Threshold sensitivity")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"c")
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,ext,"Supp11_LR_sensitivity."+ext),bbox_inches="tight",**kw)
print(perm.to_string(index=False)); print(stab.to_string(index=False)); print(dict(zip(thr,counts)))
