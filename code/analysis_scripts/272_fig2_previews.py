# -*- coding: utf-8 -*-
"""272_fig2_previews.py — six standalone previews for Figure 2."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP2=os.path.join(OUT,"step2"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig2","preview")
os.makedirs(DEST,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7.5,
 "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
 "axes.linewidth":0.6,"savefig.facecolor":"white"})
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
PHC="#B2182B"; GREY="#BBBBBB"
def save(fig,n):
    fig.savefig(os.path.join(DEST,n+".png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(DEST,n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig); print("saved",n,flush=True)

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
# a cell cycle
cc=pd.read_csv(os.path.join(OUT,"cellcycle_phase_composition.csv")).set_index("ct").reindex(["PH","GH","MH1","MH2"])
fig,ax=plt.subplots(figsize=(3.4,2.7))
cols=[PHC if c=="PH" else "#72B7B2" for c in cc.index]
ax.bar(cc.index,cc["G2M"],color=cols,width=0.62)
for i,v in enumerate(cc["G2M"]): ax.text(i,v+1,f"{v:.1f}%",ha="center",fontsize=8)
ax.set_ylabel("G2M fraction (%)"); ax.set_ylim(0,70); ax.set_title("PH is the proliferative state")
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2a_cellcycle")

# b PH markers dotplot
markers=[("PPAF2","113812814"),("Igfbp7","113811702"),("Grik2","113823802"),("SSPO","113824171"),("ccnd2","113814652")]
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
mat=[]; pct=[]
for lab,gid in markers:
    x=np.asarray(a[:,g2v[gid]].layers["lognorm"].todense()).ravel()
    mat.append([x[(a.obs.cell_type==ct).values].mean() for ct in order])
    pct.append([(x[(a.obs.cell_type==ct).values]>0).mean()*100 for ct in order])
mat=np.array(mat); pct=np.array(pct)
fig,ax=plt.subplots(figsize=(4.4,2.4))
for i in range(len(markers)):
    for j in range(len(order)):
        ax.scatter(j,i,s=max(pct[i,j]*1.3,1),c=mat[i,j],cmap="Reds",vmin=0,vmax=mat.max(),edgecolors="#666",linewidths=0.25)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=60)
ax.set_yticks(range(len(markers))); ax.set_yticklabels([m[0] for m in markers])
ax.set_title("PH progenitor markers"); ax.set_xlim(-0.6,len(order)-0.4); ax.set_ylim(len(markers)-0.5,-0.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2b_PH_markers")

# c PH substates trajectory
t=ad.read_h5ad(os.path.join(STEP2,"step2_integrated.h5ad"),backed="r")
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv")); ph.columns=["cell","group","sub","u1","u2"]
sub=ph.set_index("cell")["sub"]
d=t.obs.loc[t.obs.cell_type=="PH",["projected_pseudotime"]].copy(); d["sub"]=sub.reindex(d.index)
fig,ax=plt.subplots(figsize=(3.2,2.7)); order_s=[1,0,2,3]
for pos,s_ in enumerate(order_s):
    v=d.loc[d["sub"]==s_,"projected_pseudotime"].dropna().values
    ax.boxplot([v],positions=[pos],widths=0.5,patch_artist=True,showfliers=False,
               boxprops=dict(facecolor=PHC,alpha=0.65),medianprops=dict(color="#222",lw=0.7))
ax.set_xticks(range(4)); ax.set_xticklabels([f"s{s_}" for s_ in order_s])
ax.set_ylabel("Projected pseudotime (PH→GH)"); ax.set_title("PH substate direction")
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2c_PH_substates")

# d LyO UCell
fig,ax=plt.subplots(figsize=(3.8,2.7))
data=[a.obs.loc[a.obs.cell_type==ct,"ucell_LyO_Progenitor"].values for ct in order]
bp=ax.boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62,medianprops=dict(color="#222",lw=0.6))
for p,c in zip(bp["boxes"],order): p.set_facecolor(PHC if c=="PH" else "#72B7B2"); p.set_alpha(0.7); p.set_edgecolor("#444"); p.set_linewidth(0.4)
ax.set_ylabel("UCell score"); ax.set_title("LyO progenitor programme"); ax.tick_params(axis="x",rotation=60)
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2d_LyO_ucell")

# e overlap
ov=pd.read_csv(os.path.join(STEP3,"step3_hypergeom_overlaps.csv"))
sub=ov[ov.lyo_set=="LyO_Progenitor"].set_index("compared_to").reindex(["PH_markers","niche_markers","TS_markers","MH2_markers","MyoEpi_markers","M_markers"])
fig,ax=plt.subplots(figsize=(3.8,2.5))
y=np.arange(len(sub)); vals=sub["fold"].fillna(0).values
ax.barh(y,vals,color=[PHC if x>1 else GREY for x in vals],height=0.62)
ax.set_yticks(y); ax.set_yticklabels(["PH","niche","TS","MH2","MyoEpi","M"]); ax.invert_yaxis()
for i,(v,p,k) in enumerate(zip(vals,sub["p_hypergeom"].values,sub["overlap"].values)):
    ax.text(max(v,0.02)+2,i,f"{int(k)} genes, {v:.1f}×\nP={p:.1e}",va="center",fontsize=6)
ax.set_xlim(0,max(12,vals.max()*1.6)); ax.set_xlabel("Fold enrichment"); ax.set_title("LyO progenitor vs Hs markers")
for s in ["top","right"]: ax.spines[s].set_visible(False)
save(fig,"Fig2e_LyO_overlap")

# f eLife prohemocyte
es=pd.read_csv(os.path.join(OUT,"elife_percell_scores.csv"),index_col=0)
es=es.reindex(a.obs_names)
fig,axes=plt.subplots(1,2,figsize=(5.6,2.7),gridspec_kw={"width_ratios":[1.4,1]})
for k,(col,lab,color) in enumerate([("PH1","eLife PH1 (cycling)","#4C78A8"),("PH2","eLife PH2 (identity)","#E45756")]):
    data=[es.loc[a.obs.cell_type==ct,col].values for ct in order]
    bp=axes[0].boxplot(data,positions=np.arange(len(order))+(k-0.5)*0.3,widths=0.26,patch_artist=True,showfliers=False,medianprops=dict(color="#222",lw=0.5))
    for p in bp["boxes"]: p.set_facecolor(color); p.set_alpha(0.7); p.set_edgecolor("#444"); p.set_linewidth(0.3)
axes[0].set_xticks(range(len(order))); axes[0].set_xticklabels(order,rotation=60); axes[0].set_ylabel("eLife subtype score")
axes[0].set_title("eLife prohemocyte programmes"); axes[0].legend(handles=[Patch(color="#4C78A8",label="PH1 (cycling)"),Patch(color="#E45756",label="PH2 (identity)")],frameon=False,fontsize=5.5,loc="upper right")
ov2=pd.read_csv(os.path.join(OUT,"figures_pub","Supplementary_niche","Supp6_hypergeom.csv"))
ov2=ov2[ov2.set.isin(["eLife PH1","eLife PH2"])]
axes[1].bar(np.arange(len(ov2)),ov2.fold,color=["#4C78A8" if "PH1" in s else "#E45756" for s in ov2.set])
axes[1].set_xticks(np.arange(len(ov2))); axes[1].set_xticklabels([s.replace("eLife ","")+"\nvs "+c.replace(" markers","") for s,c in zip(ov2.set,ov2.compared_to)],rotation=45,ha="right",fontsize=6)
for i,r in enumerate(ov2.itertuples()): axes[1].text(i,r.fold+2,f"{int(r.overlap)} genes\nP={r.p:.1e}",ha="center",fontsize=5.5)
axes[1].set_ylabel("Fold enrichment"); axes[1].set_title("Overlap with Hs markers")
for ax_ in axes:
    for s in ["top","right"]: ax_.spines[s].set_visible(False)
save(fig,"Fig2f_eLife_prohemocyte")
print("done")
