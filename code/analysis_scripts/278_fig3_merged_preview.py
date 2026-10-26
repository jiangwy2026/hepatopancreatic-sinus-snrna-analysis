# -*- coding: utf-8 -*-
"""278_fig3_merged_preview.py — merged Fig3: niche landscape, LR pairs, receptors, dynamics, exclusion."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1"); STEP2=os.path.join(OUT,"step2"); STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig3_merged","preview"); os.makedirs(DEST,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
def letter(ax,s): ax.text(-0.08,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")

lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad")); T=a.obsm["X_tsne"]
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
itgb1=expr(["113830322","113805715"]); notch1=expr(["113802681"])
t2=ad.read_h5ad(os.path.join(STEP2,"step2_integrated.h5ad"),backed="r")
gdf=pd.read_csv(os.path.join(OUT,"figures_pub","Supplementary_niche","Supp7_group_LR.csv"))

fig=plt.figure(figsize=(7.2,10.2))
gs=fig.add_gridspec(4,3,height_ratios=[1.15,0.85,0.9,1.05],hspace=0.55,wspace=0.38,left=0.07,right=0.97,top=0.95,bottom=0.05)

# a circular network
ax=fig.add_subplot(gs[0,0:2]); ax.axis("off"); ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35)
srcs=["TS","MH2","MyoEpi","M"]; tgts=["PH","GH"]
ang={"TS":150,"MH2":115,"MyoEpi":80,"M":45,"PH":330,"GH":300}
cols={"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2","PH":"#E45756","GH":"#72B7B2"}
for n in srcs+tgts:
    th=np.deg2rad(ang[n]); x,y=1.02*np.cos(th),1.02*np.sin(th)
    ax.scatter(x,y,s=180 if n in tgts else 130,color=cols[n],zorder=5,edgecolors="white",linewidths=0.6)
    ax.text(1.18*np.cos(th),1.18*np.sin(th),n,ha="center",va="center",fontsize=6.5)
agg=lr.groupby(["source","target","axis"],as_index=False)["prob"].sum()
mx=agg.prob.max()
for r in agg.itertuples():
    if r.prob<0.01: continue
    a0,a1=np.deg2rad(ang[r.source]),np.deg2rad(ang[r.target])
    p0=np.array([0.98*np.cos(a0),0.98*np.sin(a0)]); p1=np.array([0.98*np.cos(a1),0.98*np.sin(a1)])
    ctrl=np.array([0.18*np.cos((a0+a1)/2),0.18*np.sin((a0+a1)/2)])
    ax.add_patch(PathPatch(Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3]),fill=False,
                           color="#B2182B" if r.axis=="ECM-integrin" else "#2166AC",lw=0.5+5*r.prob/mx,alpha=0.55,zorder=2))
ax.legend(handles=[Line2D([0],[0],color="#B2182B",lw=2,label="ECM–integrin"),Line2D([0],[0],color="#2166AC",lw=2,label="Notch")],frameon=False,fontsize=5.5,loc="lower center",bbox_to_anchor=(0.5,-0.02))
ax.set_title("Niche → hemocyte communication landscape"); letter(ax,"a")

# b key LR dot
ax=fig.add_subplot(gs[0,2])
top=lr[lr.p_perm<0.05].copy(); top["label"]=top["pair"]+"\n"+top["source"]+"→"+top["target"]
top=top.sort_values("prob",ascending=False).drop_duplicates("label").head(8).sort_values("prob")
ax.scatter(top.prob,range(len(top)),s=15+220*top.prob,c=-np.log10(top.p_perm.clip(lower=1e-300)),cmap="viridis",edgecolors="#333",linewidths=0.3,zorder=3)
for i,r in enumerate(top.itertuples()): ax.text(r.prob+0.004,i,f"P={r.p_perm:.1e}",va="center",fontsize=4.6)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.label,fontsize=4.6); ax.set_xlabel("Probability"); ax.set_xlim(0,top.prob.max()*1.35)
ax.set_title("Key ligand–receptor pairs")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"b")

# c receptor feature plots + detection
for k,(v,title) in enumerate([(itgb1,"ITGB1"),(notch1,"NOTCH1")]):
    ax=fig.add_subplot(gs[1,k])
    vmax=np.quantile(v,0.99) or 1
    sc=ax.scatter(T[:,0],T[:,1],c=v,cmap="Reds",s=1.0,linewidths=0,rasterized=True,vmin=0,vmax=vmax)
    for ct in ["PH","GH"]:
        m=(a.obs.cell_type==ct).values; ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=6,ha="center",bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title+" (receptor)")
    for s in ax.spines.values(): s.set_visible(False)
    if k==0: letter(ax,"c")
ax=fig.add_subplot(gs[1,2])
groups={"PH":["PH"],"GH":["GH"],"Niche":["TS","MH2","MyoEpi","M"],"Hepato":["F","R","B"],"Other":["PM-like"]}
pct=[(itgb1[a.obs.cell_type.isin(v).values]>0).mean()*100 for v in groups.values()]
ax.bar(range(len(groups)),pct,color="#C44E52",width=0.6)
for i,p_ in enumerate(pct): ax.text(i,p_+0.15,f"{p_:.1f}%",ha="center",fontsize=5.5)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(list(groups),rotation=40); ax.set_ylabel("ITGB1+ cells (%)"); ax.set_title("Receptor detection")
for s in ["top","right"]: ax.spines[s].set_visible(False)

# d control vs hemospasia dumbbell
ax=fig.add_subplot(gs[2,:])
piv=gdf.groupby(["source","target","group"])["prob"].sum().unstack(fill_value=0).reset_index()
piv["label"]=piv.source+"→"+piv.target; piv=piv.sort_values("Hemospasia",ascending=True)
y=np.arange(len(piv))
for i,r in piv.iterrows():
    ax.plot([r.Control,r.Hemospasia],[list(piv.index).index(i)]*2,color="#BBB",lw=1.0,zorder=1)
ax.scatter(piv.Control,range(len(piv)),s=28,c="#4C78A8",label="Control",zorder=3)
ax.scatter(piv.Hemospasia,range(len(piv)),s=28,c="#E45756",label="Hemospasia",zorder=3)
ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.label,fontsize=5.8); ax.set_xlabel("Summed ECM–integrin probability")
ax.legend(frameon=False,fontsize=6,loc="lower right"); ax.set_title("Niche signalling before vs after hemospasia")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"d")

# e trajectory
ax=fig.add_subplot(gs[3,0:2])
v=t2.obs["projected_pseudotime"].values
sc=ax.scatter(t2.obsm["X_tsne"][:,0],t2.obsm["X_tsne"][:,1],c=v,cmap="viridis",s=1.2,linewidths=0,rasterized=True)
for ct in ["PH","GH","TS","MyoEpi"]:
    m=(t2.obs.cell_type==ct).values
    ax.text(np.median(t2.obsm["X_tsne"][m,0]),np.median(t2.obsm["X_tsne"][m,1]),ct,fontsize=6,ha="center",bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
ax.set_xticks([]); ax.set_yticks([]); ax.set_title("Projection onto PH→GH trajectory")
for s in ax.spines.values(): s.set_visible(False)
plt.colorbar(sc,ax=ax,fraction=0.04,pad=0.02).set_label("Projected pseudotime",fontsize=6)
letter(ax,"e")
ax=fig.add_subplot(gs[3,2])
obs=t2.obs
data=[obs.loc[obs.cell_type=="PH","projected_pseudotime"].values,
      obs.loc[obs.cell_type=="GH","projected_pseudotime"].values,
      obs.loc[obs.compartment=="niche","projected_pseudotime"].values,
      obs.loc[obs.compartment=="hepatopancreas","projected_pseudotime"].values,
      obs.loc[obs.compartment=="other","projected_pseudotime"].values]
bp=ax.boxplot(data,tick_labels=["PH","GH","Niche","Hepato","Other"],showfliers=False,patch_artist=True,widths=0.6)
for p,c in zip(bp["boxes"],["#E45756","#72B7B2","#4C78A8","#8C510A","#999"]): p.set_facecolor(c); p.set_alpha(0.7)
ax.set_ylabel("Projected pseudotime"); ax.set_title("Compartment comparison")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.savefig(os.path.join(DEST,"Fig3_merged_preview.png"),dpi=300,bbox_inches="tight")
print("saved")
