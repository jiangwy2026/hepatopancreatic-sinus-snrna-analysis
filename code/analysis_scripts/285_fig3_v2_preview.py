# -*- coding: utf-8 -*-
"""285_fig3_v2_preview.py — readability-optimized Fig3 preview."""
import os
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig3_v2","preview"); os.makedirs(DEST,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")

lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
gdf=pd.read_csv(os.path.join(OUT,"figures_pub","Supplementary_niche","Supp7_group_LR.csv"))
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad")); T=a.obsm["X_tsne"]
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
itgb1=expr(["113830322","113805715"])

fig=plt.figure(figsize=(7.2,10.0))
gs=fig.add_gridspec(3,2,height_ratios=[1.1,0.95,1.05],hspace=0.55,wspace=0.38,left=0.08,right=0.97,top=0.95,bottom=0.05)

# A faceted heatmaps
src=["TS","MH2","MyoEpi","M"]; tgt=["PH","GH"]
sub=gs[0,0].subgridspec(1,2,wspace=0.45)
for k,axis in enumerate(["ECM-integrin","Notch"]):
    ax=fig.add_subplot(sub[0,k])
    piv=lr[lr.axis==axis].groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=src,columns=tgt,fill_value=0)
    im=ax.imshow(piv.values,cmap="Reds" if k==0 else "Blues",aspect="auto")
    ax.set_xticks([0,1]); ax.set_xticklabels(tgt); ax.set_yticks(range(4)); ax.set_yticklabels(src)
    for i in range(4):
        for j in range(2): ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",fontsize=6,color="white" if piv.values[i,j]>piv.values.max()*0.6 else "#333")
    ax.set_title(axis); ax.tick_params(length=0)
    plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02)
letter(fig.axes[0],"a")

# B lollipop
ax=fig.add_subplot(gs[0,1])
top=lr[lr.p_perm<0.05].copy().sort_values("prob",ascending=False).drop_duplicates(["pair","source","target"]).head(8)
top["label"]=top.pair+"\n"+top.source+"→"+top.target; top=top.sort_values(["source","prob"],ascending=[True,True])
y=np.arange(len(top)); cols=["#E45756" if t=="PH" else "#72B7B2" for t in top.target]
ax.hlines(y,0,top.prob,color="#DDD",lw=0.9)
ax.scatter(top.prob,y,c=cols,s=26,zorder=3)
for i,r in enumerate(top.itertuples()):
    star="***" if r.p_perm<0.001 else "**" if r.p_perm<0.01 else "*"
    ax.text(r.prob+0.004,i,star,va="center",fontsize=7)
ax.set_yticks(y); ax.set_yticklabels(top.label,fontsize=4.4); ax.set_xlabel("Communication probability"); ax.set_xlim(0,top.prob.max()*1.15)
ax.legend(handles=[Patch(color="#E45756",label="PH"),Patch(color="#72B7B2",label="GH")],frameon=False,fontsize=5.5,loc="lower right")
ax.set_title("Key ligand–receptor pairs")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"b")

# C ITGB1 violin
ax=fig.add_subplot(gs[1,0])
groups={"PH":["PH"],"GH":["GH"],"Niche":["TS","MH2","MyoEpi","M"],"Hepato":["F","R","B"],"Other":["PM-like"]}
data=[itgb1[a.obs.cell_type.isin(v).values] for v in groups.values()]
parts=ax.violinplot(data,showextrema=False,widths=0.8)
for i,b in enumerate(parts["bodies"]): b.set_facecolor("#C44E52" if i<2 else "#8DA0CB"); b.set_alpha(0.7); b.set_edgecolor("#444"); b.set_linewidth(0.4)
for i,d in enumerate(data): ax.plot([i,i],[np.median(d),np.median(d)],color="#222",lw=1.0)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(list(groups),rotation=30); ax.set_ylabel("ITGB1 (log1p CP10K)")
ax.set_title("ITGB1 receptor expression")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"c")

# D volcano
ax=fig.add_subplot(gs[1,1])
piv=gdf.groupby(["source","target","group"])["prob"].sum().unstack(fill_value=0).reset_index()
piv["label"]=piv.source+"→"+piv.target; piv["log2fc"]=np.log2((piv.Hemospasia+1e-6)/(piv.Control+1e-6))
# permutation p using cell labels
labels=a.obs.cell_type.astype(str).values; grp=a.obs.group.astype(str).replace({"blank":"Control","treated":"Hemospasia"}).values
rng=np.random.default_rng(42); ps=[]
for r in piv.itertuples():
    pairs=lr[(lr.source==r.source)&(lr.target==r.target)]
    lg=np.zeros(a.n_obs); rg=np.zeros(a.n_obs)
    for p_ in pairs.itertuples():
        lg+=expr(str(p_.ligand_genes).split(";")); rg+=expr(str(p_.receptor_genes).split(";"))
    sm=(labels==r.source); tm=(labels==r.target); obs=lg[sm&(grp=="Hemospasia")].mean()*rg[tm&(grp=="Hemospasia")].mean()-lg[sm&(grp=="Control")].mean()*rg[tm&(grp=="Control")].mean()
    null=[]
    for _ in range(500):
        perm=rng.permutation(grp)
        null.append(lg[sm&(perm=="Hemospasia")].mean()*rg[tm&(perm=="Hemospasia")].mean()-lg[sm&(perm=="Control")].mean()*rg[tm&(perm=="Control")].mean())
    null=np.array(null); p=(np.sum(np.abs(null)>=abs(obs))+1)/len(null); ps.append(p)
piv["p"]=ps
ax.scatter(piv.log2fc,-np.log10(piv.p.clip(lower=1e-4)),s=34,c=["#C44E52" if p<0.1 else "#999" for p in piv.p],edgecolors="#333",linewidths=0.3,zorder=3)
for r in piv.itertuples():
    if r.p<0.1: ax.annotate(r.label,(r.log2fc,-np.log10(max(r.p,1e-4))),textcoords="offset points",xytext=(4,3),fontsize=5.2)
ax.axhline(-np.log10(0.1),color="#999",ls="--",lw=0.5); ax.axvline(0,color="#999",ls="--",lw=0.5)
ax.set_xlabel("log2 FC (Hemospasia / Control)"); ax.set_ylabel("−log10 P"); ax.set_title("Dynamic niche signalling")
for s in ["top","right"]: ax.spines[s].set_visible(False); letter(ax,"d")

# E ridge + box
ax=fig.add_subplot(gs[2,0])
obs=a.obs.copy(); obs["grp2"]=obs.cell_type.astype(str).replace({"F":"Hepato","R":"Hepato","B":"Hepato","PM-like":"Other"})
cats=["PH","GH","MyoEpi","M","TS","MH2","Hepato","Other"]
xs=np.linspace(0,0.55,220)
for i,c in enumerate(cats[::-1]):
    v=obs.loc[obs.grp2==c,"projected_pseudotime"].dropna().values
    if len(v)<10: continue
    y=gaussian_kde(v)(xs); y=y/y.max()*0.8
    ax.fill_between(xs,i,i+y,color="#B2182B" if c=="PH" else ("#72B7B2" if c=="GH" else "#8DA0CB"),alpha=0.75,lw=0)
    ax.text(0.002,i+0.12,c,ha="left",va="center",fontsize=5.8)
ax.set_yticks([]); ax.set_xlim(0,0.5); ax.set_xlabel("Projected pseudotime"); ax.set_title("Pseudotime distribution by state")
for s in ["top","right","left"]: ax.spines[s].set_visible(False); letter(ax,"e")

ax=fig.add_subplot(gs[2,1])
data=[obs.loc[obs.grp2==c,"projected_pseudotime"].dropna().values for c in cats]
bp=ax.boxplot(data,tick_labels=cats,showfliers=False,patch_artist=True,widths=0.6)
for p,c in zip(bp["boxes"],["#B2182B" if x=="PH" else ("#72B7B2" if x=="GH" else "#8DA0CB") for x in cats]): p.set_facecolor(c); p.set_alpha(0.7)
ax.set_ylabel("Projected pseudotime"); ax.set_title("Compartment comparison"); ax.tick_params(axis="x",rotation=60)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.savefig(os.path.join(DEST,"Fig3_v2_preview.png"),dpi=300,bbox_inches="tight")
print("saved")
