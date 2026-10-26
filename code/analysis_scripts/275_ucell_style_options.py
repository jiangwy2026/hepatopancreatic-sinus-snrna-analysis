# -*- coding: utf-8 -*-
"""275_ucell_style_options.py — violin / ridge / lollipop options for UCell."""
import os, collections
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig2","preview")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7.5,
 "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
 "axes.linewidth":0.6,"savefig.facecolor":"white"})
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
PHC="#B2182B"; GREY="#8DA0CB"

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
vn=set(map(str,a.var_names)); sym=collections.defaultdict(list)
for i,s in enumerate(a.var["sym"].astype(str)): sym[s.lower()].append(a.var_names[i])
def map_set(df,n=50):
    d=df.sort_values("Log2FC",ascending=False).head(n); out=[]
    for x in d["Gene Name"].astype(str):
        v=x if x in vn else (sym[x.lower()][0] if len(sym.get(x.lower(),[]))==1 else None)
        if v and v not in out: out.append(v)
    return out
p1=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp3-v2.csv")),50)
p2=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp6-v2.csv")),50)
genes=list(dict.fromkeys(p1+p2)); idx=[list(a.var_names).index(g) for g in genes]
N=a.n_vars; k=len(idx); score=np.zeros(a.n_obs,dtype=np.float32)
for st in range(0,a.n_obs,400):
    en=min(a.n_obs,st+400); X=np.asarray(a[st:en].X.todense(),dtype=np.float32)
    rank=np.argsort(-X,axis=1,kind="stable"); r=np.empty_like(rank); r[np.arange(en-st)[:,None],rank]=np.arange(1,N+1)
    S=r[:,idx].sum(1); score[st:en]=1-(S-k*(k+1)/2)/(k*(N-k))
a.obs["blood_prog"]=score

# violin
fig,ax=plt.subplots(figsize=(5.2,2.8))
data=[a.obs.loc[a.obs.cell_type==ct,"blood_prog"].values for ct in order]
parts=ax.violinplot(data,showextrema=False,widths=0.8)
for i,b in enumerate(parts["bodies"]):
    b.set_facecolor(PHC if order[i]=="PH" else GREY); b.set_alpha(0.75); b.set_edgecolor("#444"); b.set_linewidth(0.4)
for i,v in enumerate(data):
    ax.plot([i,i],[np.median(v),np.median(v)],color="#222",lw=1.0)
    ax.text(i,ax.get_ylim()[1] if False else 0.668,f"n={len(v)}",ha="center",fontsize=5)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=60); ax.set_ylabel("UCell score")
ax.set_title("Blood-cell prohemocyte programme (violin)"); ax.set_ylim(0.57,0.67)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(DEST,"Fig2f_style_violin.png"),dpi=300,bbox_inches="tight"); plt.close(fig)

# ridge
fig,ax=plt.subplots(figsize=(5.0,3.0))
xs=np.linspace(0.56,0.68,200)
for i,ct in enumerate(order[::-1]):
    v=a.obs.loc[a.obs.cell_type==ct,"blood_prog"].values
    y=gaussian_kde(v)(xs); y=y/y.max()*0.85
    ax.fill_between(xs,i,i+y,color=PHC if ct=="PH" else GREY,alpha=0.8,lw=0)
    ax.plot(xs,i+y,color="#333",lw=0.4)
    ax.text(0.558,i+0.15,ct,ha="right",va="center",fontsize=6)
ax.set_yticks([]); ax.set_xlim(0.555,0.68); ax.set_xlabel("UCell score")
ax.set_title("Blood-cell prohemocyte programme (ridge)")
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(DEST,"Fig2f_style_ridge.png"),dpi=300,bbox_inches="tight"); plt.close(fig)

# lollipop
fig,ax=plt.subplots(figsize=(4.4,2.7))
means=[np.mean(v) for v in data]; pct=[np.mean(v>np.median(score))*100 for v in data]
ax.hlines(range(len(order)),0.57,means,color="#CCC",lw=0.8)
sc=ax.scatter(means,range(len(order)),s=[20+p*2 for p in pct],c=[PHC if ct=="PH" else GREY for ct in order],edgecolors="#333",linewidths=0.4,zorder=3)
for i,(m,p) in enumerate(zip(means,pct)): ax.text(m+0.002,i,f"{m:.3f} ({p:.0f}%)",va="center",fontsize=5.5)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
ax.set_xlim(0.57,0.65); ax.set_xlabel("Mean UCell score"); ax.set_title("Blood-cell prohemocyte programme (lollipop)")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(DEST,"Fig2f_style_lollipop.png"),dpi=300,bbox_inches="tight"); plt.close(fig)
print("done")
