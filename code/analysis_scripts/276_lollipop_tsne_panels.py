# -*- coding: utf-8 -*-
"""276_lollipop_tsne_panels.py — t-SNE + lollipop for LyO and blood-cell UCell."""
import os, collections
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig2","preview")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7.5,
 "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
 "axes.linewidth":0.6,"savefig.facecolor":"white","svg.fonttype":"none"})
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
PHC="#B2182B"; GREY="#8DA0CB"

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
T=a.obsm["X_tsne"]
# combined blood-cell programme
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
N=a.n_vars; k=len(idx); blood=np.zeros(a.n_obs,dtype=np.float32)
for st in range(0,a.n_obs,400):
    en=min(a.n_obs,st+400); X=np.asarray(a[st:en].X.todense(),dtype=np.float32)
    rank=np.argsort(-X,axis=1,kind="stable"); r=np.empty_like(rank); r[np.arange(en-st)[:,None],rank]=np.arange(1,N+1)
    S=r[:,idx].sum(1); blood[st:en]=1-(S-k*(k+1)/2)/(k*(N-k))

def panel(score,title,fname):
    thr=np.median(score)
    fig=plt.figure(figsize=(7.2,3.0)); gs=fig.add_gridspec(1,2,width_ratios=[1.15,1],wspace=0.35,left=0.07,right=0.97,top=0.88,bottom=0.22)
    ax=fig.add_subplot(gs[0])
    sca=ax.scatter(T[:,0],T[:,1],c=score,cmap="viridis",s=1.4,linewidths=0,rasterized=True,
                   vmin=np.quantile(score,0.02),vmax=np.quantile(score,0.98))
    for ct in ["PH","GH","TS"]:
        m=(a.obs.cell_type==ct).values
        ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=6.5,ha="center",
                bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.7))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_title(title+" (t-SNE)")
    for s in ax.spines.values(): s.set_visible(False)
    cb=plt.colorbar(sca,ax=ax,fraction=0.045,pad=0.02); cb.set_label("UCell score",fontsize=7); cb.ax.tick_params(labelsize=6)
    ax2=fig.add_subplot(gs[1])
    means=[score[(a.obs.cell_type==ct).values].mean() for ct in order]
    pct=[(score[(a.obs.cell_type==ct).values]>thr).mean()*100 for ct in order]
    ax2.hlines(range(len(order)),min(means)*0.99,means,color="#DDD",lw=0.8)
    ax2.scatter(means,range(len(order)),s=[18+p*1.8 for p in pct],c=[PHC if ct=="PH" else GREY for ct in order],edgecolors="#333",linewidths=0.4,zorder=3)
    for i,(m,p) in enumerate(zip(means,pct)): ax2.text(m+0.0025,i,f"{m:.3f} ({p:.0f}%)",va="center",fontsize=5.5)
    ax2.set_yticks(range(len(order))); ax2.set_yticklabels(order); ax2.invert_yaxis()
    ax2.set_xlim(min(means)*0.985,max(means)*1.03); ax2.set_xlabel("Mean UCell score")
    ax2.set_title(title+" (summary)")
    for s in ["top","right"]: ax2.spines[s].set_visible(False)
    fig.savefig(os.path.join(DEST,fname+".png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(DEST,fname+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,fname+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig); print("saved",fname,flush=True)

panel(a.obs["ucell_LyO_Progenitor"].values,"LyO progenitor programme","Fig2d_LyO_ucell_lollipop_tsne")
panel(blood,"Blood-cell prohemocyte programme","Fig2f_blood_prog_lollipop_tsne")
print("done")
