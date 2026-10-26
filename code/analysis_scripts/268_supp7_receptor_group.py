# -*- coding: utf-8 -*-
"""268_supp7_receptor_group.py — receptor expression + Control/Hemospasia LR comparison."""
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
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
T=a.obsm["X_tsne"]; g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def mean_expr(ids):
    idx=[g2v[i] for i in ids if i in g2v]
    return np.asarray(a[:,idx].layers["lognorm"].todense()).mean(1) if idx else np.zeros(a.n_obs)
expr={"ITGB1":mean_expr(["113830322","113805715"]),
      "ITGAV":mean_expr(["113800203","113821319"]),
      "ITGA4":mean_expr(["113827435"]),
      "NOTCH1":mean_expr(["113802681"])}

fig=plt.figure(figsize=(7.4,4.8)); gs=fig.add_gridspec(1,3,width_ratios=[1.2,1.0,1.25],wspace=0.48,left=0.07,right=0.97,top=0.9,bottom=0.2)
sub=gs[0].subgridspec(2,2,hspace=0.25,wspace=0.25)
for k,(name,v) in enumerate(expr.items()):
    ax=fig.add_subplot(sub[k//2,k%2])
    vmax=np.quantile(v,0.99); vmax=vmax if vmax>0 else 1
    ax.scatter(T[:,0],T[:,1],c=v,cmap="Reds",s=1.0,linewidths=0,rasterized=True,vmin=0,vmax=vmax)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(name,fontsize=6.5)
    for sp in ax.spines.values(): sp.set_visible(False)
    if k==0: letter(ax,"a")

ax=fig.add_subplot(gs[1])
groups={"PH":["PH"],"GH":["GH"],"Niche":["TS","MH2","MyoEpi","M"],"Hepato":["F","R","B"],"Other":["PM-like"]}
pct=[]; mpos=[]
for v in groups.values():
    x=expr["ITGB1"][a.obs.cell_type.isin(v).values]
    pct.append((x>0).mean()*100); mpos.append(x[x>0].mean() if (x>0).any() else 0)
ax.bar(range(len(groups)),pct,color="#C44E52",width=0.62)
for i,(p_,m_) in enumerate(zip(pct,mpos)): ax.text(i,p_+0.2,f"{p_:.1f}%\n(mean+ {m_:.2f})",ha="center",fontsize=4.8)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(list(groups),rotation=35,fontsize=5.5)
ax.set_ylabel("ITGB1-positive cells (%)"); ax.set_title("ITGB1 detection in receiver states")
ax.set_ylim(0,max(pct)*1.3)
for sp in ["top","right"]: ax.spines[sp].set_visible(False)
letter(ax,"b")

lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
lr=lr[lr.axis=="ECM-integrin"].copy()
gcol=a.obs["group"].astype(str).replace({"blank":"Control","treated":"Hemospasia"})
rows=[]
for r in lr.itertuples():
    lg=mean_expr(str(r.ligand_genes).split(";")) if pd.notna(r.ligand_genes) else np.zeros(a.n_obs)
    rg=mean_expr(str(r.receptor_genes).split(";")) if pd.notna(r.receptor_genes) else np.zeros(a.n_obs)
    for grp in ["Control","Hemospasia"]:
        m=(gcol==grp).values
        rows.append({"pair":r.pair,"source":r.source,"target":r.target,"group":grp,
                     "prob":float(lg[m & a.obs.cell_type.eq(r.source).values].mean()*rg[m & a.obs.cell_type.eq(r.target).values].mean())})
gdf=pd.DataFrame(rows); gdf.to_csv(os.path.join(DEST,"Supp7_group_LR.csv"),index=False,encoding="utf-8-sig")
ax=fig.add_subplot(gs[2])
piv={grp:gdf[gdf.group==grp].groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=["TS","MH2","MyoEpi","M"],columns=["PH","GH"],fill_value=0) for grp in ["Control","Hemospasia"]}
mat=np.concatenate([piv["Control"].values,piv["Hemospasia"].values],axis=1)
im=ax.imshow(mat,cmap="Reds",aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["PH\nC","GH\nC","PH\nH","GH\nH"],fontsize=5.2)
ax.set_yticks(range(4)); ax.set_yticklabels(["TS","MH2","MyoEpi","M"])
for i in range(4):
    for j in range(4): ax.text(j,i,f"{mat[i,j]:.2f}",ha="center",va="center",fontsize=5,color="white" if mat[i,j]>mat.max()*0.6 else "#333")
ax.set_title("ECM-integrin: Control vs Hemospasia"); plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02)
letter(ax,"c")
fig.savefig(os.path.join(DEST,"png","Supp7_receptor_group.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(DEST,"svg","Supp7_receptor_group.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"tiff","Supp7_receptor_group.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
print("done")
