# -*- coding: utf-8 -*-
"""273_fig2_revised_previews.py — revised 2a, 2e, 2f."""
import os, collections
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import hypergeom
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
 "axes.linewidth":0.6,"savefig.facecolor":"white","svg.fonttype":"none"})
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
PHC="#B2182B"; GREY="#BBBBBB"
def save(fig,n):
    fig.savefig(os.path.join(DEST,n+".png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(DEST,n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig); print("saved",n,flush=True)

# 2a cell cycle without MH1
cc=pd.read_csv(os.path.join(OUT,"cellcycle_phase_composition.csv")).set_index("ct")
dyn=pd.read_csv(os.path.join(OUT,"Fig5_dynamic_cellcycle_by_group.csv"))
pm=dyn[dyn.population=="PM-like"]
pm_g2m=(pm.n*pm.pct_G2M).sum()/pm.n.sum()
states=["PH","GH","MH2","PM-like"]; vals=[cc.loc["PH","G2M"],cc.loc["GH","G2M"],cc.loc["MH2","G2M"],pm_g2m]
fig,ax=plt.subplots(figsize=(3.2,2.7))
ax.bar(states,vals,color=[PHC,"#72B7B2","#F58518","#79706E"],width=0.62)
for i,v in enumerate(vals): ax.text(i,v+1,f"{v:.1f}%",ha="center",fontsize=8)
ax.set_ylabel("G2M fraction (%)"); ax.set_ylim(0,70); ax.set_title("PH is the proliferative state")
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2a_cellcycle_v2")

# combined blood-cell prohemocyte programme
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
gene_set=list(dict.fromkeys(p1+p2)); idx=[list(a.var_names).index(g) for g in gene_set]
N=a.n_vars; k=len(idx); score=np.zeros(a.n_obs,dtype=np.float32)
for st in range(0,a.n_obs,400):
    en=min(a.n_obs,st+400); X=np.asarray(a[st:en].X.todense(),dtype=np.float32)
    rank=np.argsort(-X,axis=1,kind="stable"); r=np.empty_like(rank); r[np.arange(en-st)[:,None],rank]=np.arange(1,N+1)
    S=r[:,idx].sum(1); score[st:en]=1-(S-k*(k+1)/2)/(k*(N-k))
a.obs["blood_prog"]=score
print("combined genes",k)
# 2e t-SNE
T=a.obsm["X_tsne"]
fig,ax=plt.subplots(figsize=(3.6,3.1))
sc=ax.scatter(T[:,0],T[:,1],c=score,cmap="viridis",s=1.4,linewidths=0,rasterized=True,
              vmin=np.quantile(score,0.02),vmax=np.quantile(score,0.98))
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
ax.set_title("Blood-cell prohemocyte programme")
plt.colorbar(sc,ax=ax,fraction=0.045,pad=0.02).set_label("UCell score",fontsize=7)
for s in ax.spines.values(): s.set_visible(False)
save(fig,"Fig2e_blood_program_tsne")
# 2f boxplot + overlap
fig,axes=plt.subplots(1,2,figsize=(5.8,2.8),gridspec_kw={"width_ratios":[1.35,1]})
data=[a.obs.loc[a.obs.cell_type==ct,"blood_prog"].values for ct in order]
bp=axes[0].boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62,medianprops=dict(color="#222",lw=0.6))
for p,c in zip(bp["boxes"],order): p.set_facecolor(PHC if c=="PH" else "#72B7B2"); p.set_alpha(0.7); p.set_edgecolor("#444"); p.set_linewidth(0.4)
axes[0].set_ylabel("UCell score"); axes[0].set_title("Blood-cell prohemocyte programme"); axes[0].tick_params(axis="x",rotation=60)
mk=pd.read_csv(os.path.join(OUT,"markers_10groups.csv"))
def top(grp,n=50):
    g=mk[mk.grp==grp].sort_values("scores",ascending=False).head(n)
    return set(x for x in g["names"].astype(str) if x in vn)
ph=top("PH"); niche=set().union(*[top(g) for g in ["TS","MH2","MyoEpi","M"]])
rows=[]
for lab,B in [("PH",ph),("niche",niche)]:
    kk=len(set(gene_set)&B); M=len(vn); n=len(gene_set); N2=len(B)
    rows.append((lab,kk,kk/(n*N2/M),hypergeom.sf(kk-1,M,n,N2) if kk else 1.0))
axes[1].bar([r[0] for r in rows],[r[2] for r in rows],color=[PHC,GREY],width=0.55)
mx=max(r[2] for r in rows)
axes[1].set_ylim(0,mx*1.35)
for i,r in enumerate(rows):
    if r[2]>mx*0.5: axes[1].text(i,r[2]*0.88,f"{r[1]} genes\nP={r[3]:.1e}",ha="center",va="top",fontsize=6,color="white")
    else: axes[1].text(i,1,f"{r[1]} genes\nP={r[3]:.1e}",ha="center",va="bottom",fontsize=6,color="#333")
axes[1].set_ylabel("Fold enrichment"); axes[1].set_title("Overlap with Hs markers")
for ax_ in axes:
    for s in ["top","right"]: ax_.spines[s].set_visible(False)
save(fig,"Fig2f_blood_program_combined")
print("done")
