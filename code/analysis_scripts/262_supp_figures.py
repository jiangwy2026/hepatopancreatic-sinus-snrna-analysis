# -*- coding: utf-8 -*-
"""262_supp_figures.py — supplementary figures for the niche / cross-organ part.
All figures: Times New Roman, SVG + 1200-dpi TIFF + PNG preview.
"""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import ListedColormap, BoundaryNorm

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
for s in ["svg","tiff","png"]: os.makedirs(os.path.join(DEST,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})

def save(fig,name):
    fig.savefig(os.path.join(DEST,"svg",name+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",name+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",name+".png"),dpi=200,bbox_inches="tight")
    plt.close(fig); print("saved",name,flush=True)

def letter(ax,s): ax.text(-0.10,1.04,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]

# S1 niche marker dotplot
genes=[("113824171","VWF/hemocytin"),("113823118","THBS4"),("113830322","ITGB1 (a)"),
       ("113805715","ITGB1 (b)"),("113821319","ITGAV"),("113827435","ITGA4"),
       ("113821858","JAG1"),("113802681","NOTCH1"),("113822052","Astakine")]
rows=[]
for gid,lab in genes:
    if gid not in g2v: continue
    x=np.asarray(a[:,g2v[gid]].layers["lognorm"].todense()).ravel()
    for ct in order:
        m=(a.obs.cell_type==ct).values
        rows.append({"gene":lab,"cell_type":ct,"mean":float(x[m].mean()),"pct":float((x[m]>0).mean()*100)})
d=pd.DataFrame(rows)
fig,ax=plt.subplots(figsize=(5.4,3.2))
for i,lab in enumerate([g[1] for g in genes if g[0] in g2v]):
    for j,ct in enumerate(order):
        r=d[(d.gene==lab)&(d.cell_type==ct)].iloc[0]
        ax.scatter(j,i,s=max(r.pct*1.4,1),c=r["mean"],cmap="Reds",vmin=0,vmax=d["mean"].max(),edgecolors="#777",linewidths=0.3)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=60)
ax.set_yticks(range(len([g for g in genes if g[0] in g2v]))); ax.set_yticklabels([g[1] for g in genes if g[0] in g2v])
ax.set_title("Niche / support factors")
for v in [0,25,50,75]:
    ax.scatter([],[],s=max(v*1.4,1),c="#BBBBBB",edgecolors="#777",linewidths=0.3,label=f"{v}%")
ax.legend(title="% expressing",frameon=False,fontsize=5.5,title_fontsize=5.5,loc="center left",bbox_to_anchor=(1.01,0.5))
ax.set_xlim(-0.6,len(order)-0.4); ax.set_ylim(len(d.gene.unique())-0.5,-0.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"a"); fig.tight_layout(); save(fig,"Supp1_niche_marker_dotplot")

# S2 pair-level LR dotplot
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
pairs=lr.groupby("pair")["prob"].max().sort_values(ascending=False).head(7).index.tolist()
combos=[f"{s}->{t}" for s in ["TS","MH2","MyoEpi","M"] for t in ["PH","GH"]]
fig,ax=plt.subplots(figsize=(5.2,2.9))
for i,p in enumerate(pairs):
    for j,c in enumerate(combos):
        s,t=c.split("->"); r=lr[(lr.pair==p)&(lr.source==s)&(lr.target==t)]
        if len(r)==0: continue
        r=r.iloc[0]
        ax.scatter(j,i,s=8+900*r.prob,c=-np.log10(max(r.p_perm,1e-300)),cmap="viridis",vmin=0,vmax=4,edgecolors="#666",linewidths=0.2)
ax.set_xticks(range(len(combos))); ax.set_xticklabels(combos,rotation=60,fontsize=5.8)
ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs,fontsize=6)
ax.set_title("Pair-level niche -> hemocyte communication")
sm=plt.cm.ScalarMappable(cmap="viridis",norm=BoundaryNorm([0,1,2,3,4],256)); sm.set_array([])
cb=plt.colorbar(sm,ax=ax,fraction=0.04,pad=0.02); cb.set_label("-log10(P)",fontsize=6); cb.ax.tick_params(labelsize=5.5)
for v in [0.05,0.1,0.2]:
    ax.scatter([],[],s=8+900*v,c="#BBBBBB",edgecolors="#666",linewidths=0.2,label=f"{v:.2f}")
ax.legend(title="Probability",frameon=False,fontsize=5.5,title_fontsize=5.5,loc="center left",bbox_to_anchor=(1.18,0.5))
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"a"); fig.tight_layout(); save(fig,"Supp2_LR_pair_dotplot")

# S3 SonicParanoid presence/absence
sv=pd.read_csv(os.path.join(STEP1,"sonicparanoid_final_verdict.csv"))
neg=["CXCL12","KITLG","KIT","LYVE1","PTN","CXCR4","TEK","TIE1","CSF1R","CD34","CXCL8","CCL2","CDH5","PECAM1","FZD1","FZD7","NOTCH4","IL6","LIF","TGFB1"]
pos=["JAG1","JAG2","DLL1","DLL4","NOTCH1","NOTCH2","NOTCH3","WNT5A","FZD4","FZD5","THBS4","LAMA1","ITGAV","ITGB1","ITGA4","MRC1","VWF","POSTN","DAB2","PROM1","ANGPT1","ANGPT2","BMP4","VEGFA","FLT1","KDR","DKK1","SFRP1"]
gl=neg+pos; M=[]
for g in gl:
    r=sv[sv.human_gene==g]
    M.append(1 if (len(r) and bool(r.iloc[0].in_orthogroup)) else 0)
M=np.array(M).reshape(-1,1)
fig,ax=plt.subplots(figsize=(1.3,6.2))
ax.imshow(M,aspect="auto",cmap=ListedColormap(["#E8E8E8","#2A9D8F"]),vmin=0,vmax=1)
ax.set_xticks([0]); ax.set_xticklabels(["Orthogroup"],rotation=90,fontsize=6)
ax.set_yticks(range(len(gl))); ax.set_yticklabels(gl,fontsize=5.8)
for i,g in enumerate(gl):
    r=sv[sv.human_gene==g]
    if len(r) and bool(r.iloc[0].in_orthogroup):
        ax.text(0,i,str(int(r.iloc[0].group_id)),ha="center",va="center",fontsize=4.5,color="white")
ax.axhline(len(neg)-0.5,color="black",lw=0.6)
letter(ax,"a"); fig.tight_layout(); save(fig,"Supp3_orthogroup_presence")

# S4 per-type QC + QC-standardized composition
qc=pd.read_csv(os.path.join(OUT,"pertype_qc_stats.csv")).set_index("ct")
qs=qc.loc[[c for c in ["PH","GH","MH1","MH2","F","R","B","MyoEpi","M"] if c in qc.index],["mt_med","ng_med","umi_med"]]
Z=(qs-qs.mean())/(qs.std()+1e-9)
fig,axes=plt.subplots(1,2,figsize=(7.8,3.0))
im=axes[0].imshow(Z.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
axes[0].set_xticks([0,1,2]); axes[0].set_xticklabels(["%MT","Genes","UMI"]); axes[0].set_yticks(range(len(Z))); axes[0].set_yticklabels(Z.index)
axes[0].set_title("Per-cell-type QC medians (z-score)")
plt.colorbar(im,ax=axes[0],fraction=0.04,pad=0.02)
comp=pd.read_csv(os.path.join(OUT,"Fig5_qc_standardized_composition.csv")).set_index("cluster")
comp=comp.loc[[c for c in ["PH","GH","TS","MH2","F","R","B","MyoEpi","M","PM-like"] if c in comp.index]]
y=np.arange(len(comp))
axes[1].scatter(comp.log2fc_raw,y,s=14,c="#4C78A8",label="Raw")
axes[1].scatter(comp.log2fc_qc_std,y,s=14,c="#E45756",label="QC-standardized")
for i in range(len(comp)): axes[1].plot([comp.log2fc_raw.iloc[i],comp.log2fc_qc_std.iloc[i]],[i,i],color="#999",lw=0.5)
axes[1].axvline(0,color="black",lw=0.5,ls="--"); axes[1].set_yticks(y); axes[1].set_yticklabels(comp.index)
axes[1].set_xlabel("log2 fold change (Hemospasia / Control)"); axes[1].set_title("Composition: raw vs QC-standardized")
axes[1].legend(frameon=False,fontsize=6,loc="center left",bbox_to_anchor=(1.02,0.5))
for s in ["top","right"]: axes[1].spines[s].set_visible(False)
letter(axes[0],"a"); letter(axes[1],"b"); fig.tight_layout(); save(fig,"Supp4_QC_composition")

# S5 PH substates t-SNE + markers (computed directly from the PH subset)
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv"))
ph.columns=["cell","group","subcluster","umap1","umap2"]
sub_ids=ph.set_index("cell")["subcluster"]
mask=np.isin(a.obs_names, sub_ids.index)
A=a[mask].copy()
sub=sub_ids.reindex(A.obs_names).values
import scanpy as sc
sc.pp.highly_variable_genes(A, n_top_genes=1000)
A=A[:,A.var.highly_variable.values].copy()
sc.pp.scale(A,max_value=10); sc.tl.pca(A,n_comps=20,random_state=0)
sc.pp.neighbors(A,n_neighbors=20,random_state=0); sc.tl.tsne(A,n_pcs=20,perplexity=30,random_state=0)
T=A.obsm["X_tsne"]
fig,axes=plt.subplots(1,2,figsize=(7.6,3.1),gridspec_kw={"width_ratios":[1,1.35]})
cols={0:"#4C78A8",1:"#F58518",2:"#54A24B",3:"#B279A2"}
for s_ in sorted(set(sub)):
    m=sub==s_
    axes[0].scatter(T[m,0],T[m,1],s=2.2,c=cols.get(s_,"#999"),label=f"PH_s{s_}",rasterized=True)
axes[0].set_xticks([]); axes[0].set_yticks([]); axes[0].set_title("PH substates (t-SNE)")
axes[0].legend(frameon=False,fontsize=5.5,markerscale=2.2,loc="best")
for sp in axes[0].spines.values(): sp.set_visible(False)
# direct marker computation on PH cells
B=a[mask].copy(); B.X=B.layers["lognorm"].copy()
Xb=np.asarray(B.X.todense(),dtype=np.float32)
mt=B.var["mt"].astype(bool).values if "mt" in B.var else np.zeros(B.n_vars,bool)
syms=B.var["sym"].astype(str).values
rows=[]
for s_ in sorted(set(sub)):
    m=sub==s_
    mean_in=Xb[m].mean(0); mean_out=Xb[~m].mean(0); pct=(Xb[m]>0).mean(0)
    score=mean_in-mean_out
    ok=(~mt)&(pct>0.10)
    idx=np.argsort(np.where(ok,score,-1))[::-1][:3]
    for i in idx:
        rows.append({"subcluster":s_,"gid":B.var_names[i],"label":(syms[i] if not str(syms[i]).startswith("ncbi_") else B.var_names[i]),
                     "mean":float(mean_in[i]),"pct":float(pct[i]*100)})
top=pd.DataFrame(rows)
for i,row in enumerate(top.itertuples()):
    for j,s_ in enumerate(sorted(set(sub))):
        m=sub==s_
        x=Xb[m,:][:,list(B.var_names).index(row.gid)]
        axes[1].scatter(j,i,s=max((x>0).mean()*140,1),c=x.mean(),cmap="Reds",vmin=0,vmax=max(top["mean"].max(),1e-6),edgecolors="#777",linewidths=0.2)
axes[1].set_xticks(range(len(sorted(set(sub))))); axes[1].set_xticklabels([f"PH_s{s_}" for s_ in sorted(set(sub))])
axes[1].set_yticks(range(len(top))); axes[1].set_yticklabels([f"{r.label} (s{int(r.subcluster)})" for r in top.itertuples()],fontsize=5.2)
axes[1].set_title("Top markers per PH substate")
axes[1].set_ylim(len(top)-0.5,-0.5)
for sp in ["top","right"]: axes[1].spines[sp].set_visible(False)
letter(axes[0],"a"); letter(axes[1],"b"); fig.tight_layout(); save(fig,"Supp5_PH_substates")

print("Supplementary figures done ->",DEST)
