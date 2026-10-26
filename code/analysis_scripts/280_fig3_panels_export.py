# -*- coding: utf-8 -*-
"""280_fig3_panels_export.py — standalone Fig3 panels, SVG + 1200 dpi TIFF."""
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
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig3_merged")
for s in ["SVG","TIFF","PNG"]: os.makedirs(os.path.join(DEST,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
def letter(ax,s): ax.text(-0.09,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")
def save(fig,n):
    fig.savefig(os.path.join(DEST,"SVG",n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"TIFF",n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"PNG",n+".png"),dpi=300,bbox_inches="tight"); plt.close(fig); print("saved",n,flush=True)

lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad")); T=a.obsm["X_tsne"]
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
itgb1=expr(["113830322","113805715"]); notch1=expr(["113802681"])
t2=ad.read_h5ad(os.path.join(STEP2,"step2_integrated.h5ad"),backed="r")
gdf=pd.read_csv(os.path.join(OUT,"figures_pub","Supplementary_niche","Supp7_group_LR.csv"))

# a network
fig,ax=plt.subplots(figsize=(3.9,3.7)); ax.axis("off"); ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35)
srcs=["TS","MH2","MyoEpi","M"]; tgts=["PH","GH"]; ang={"TS":150,"MH2":115,"MyoEpi":80,"M":45,"PH":330,"GH":300}
cols={"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2","PH":"#E45756","GH":"#72B7B2"}
for n in srcs+tgts:
    th=np.deg2rad(ang[n]); ax.scatter(1.02*np.cos(th),1.02*np.sin(th),s=170 if n in tgts else 120,color=cols[n],zorder=5,edgecolors="white",linewidths=0.6)
    ax.text(1.18*np.cos(th),1.18*np.sin(th),n,ha="center",va="center",fontsize=6.5)
agg=lr.groupby(["source","target","axis"],as_index=False)["prob"].sum(); mx=agg.prob.max()
for r in agg.itertuples():
    if r.prob<0.01: continue
    a0,a1=np.deg2rad(ang[r.source]),np.deg2rad(ang[r.target]); p0=np.array([0.98*np.cos(a0),0.98*np.sin(a0)]); p1=np.array([0.98*np.cos(a1),0.98*np.sin(a1)])
    ctrl=np.array([0.18*np.cos((a0+a1)/2),0.18*np.sin((a0+a1)/2)])
    ax.add_patch(PathPatch(Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3]),fill=False,color="#B2182B" if r.axis=="ECM-integrin" else "#2166AC",lw=0.5+5*r.prob/mx,alpha=0.55,zorder=2))
ax.legend(handles=[Line2D([0],[0],color="#B2182B",lw=2,label="ECM–integrin"),Line2D([0],[0],color="#2166AC",lw=2,label="Notch")],frameon=False,fontsize=5.5,loc="lower center",bbox_to_anchor=(0.5,-0.02))
ax.set_title("Niche → hemocyte communication"); letter(ax,"a"); fig.tight_layout(); save(fig,"a_communication_network")

# b LR dot
top=lr[lr.p_perm<0.05].copy(); top["label"]=top["pair"]+"  "+top["source"]+"→"+top["target"]
top=top.sort_values("prob",ascending=False).drop_duplicates("label").head(8).sort_values("prob")
fig,ax=plt.subplots(figsize=(3.8,3.0))
ax.scatter(top.prob,range(len(top)),s=15+220*top.prob,c=-np.log10(top.p_perm.clip(lower=1e-300)),cmap="viridis",edgecolors="#333",linewidths=0.3,zorder=3)
for i,r in enumerate(top.itertuples()): ax.text(r.prob+0.004,i,f"P={r.p_perm:.1e}",va="center",fontsize=4.8)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.label,fontsize=5); ax.set_xlabel("Probability"); ax.set_xlim(0,top.prob.max()*1.35)
ax.set_title("Key ligand–receptor pairs")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"b"); fig.tight_layout(); save(fig,"b_LR_pairs")

# c1/c2 features
for v,name,lab in [(itgb1,"c1_ITGB1_feature","ITGB1"),(notch1,"c2_NOTCH1_feature","NOTCH1")]:
    fig,ax=plt.subplots(figsize=(3.4,3.0)); vmax=np.quantile(v,0.99) or 1
    sc=ax.scatter(T[:,0],T[:,1],c=v,cmap="Reds",s=1.0,linewidths=0,rasterized=True,vmin=0,vmax=vmax)
    for ct in ["PH","GH"]:
        m=(a.obs.cell_type==ct).values; ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=6,ha="center",bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(lab+" (receptor)")
    for s in ax.spines.values(): s.set_visible(False)
    plt.colorbar(sc,ax=ax,fraction=0.045,pad=0.02).set_label("log1p CP10K",fontsize=6)
    letter(ax,"c"); fig.tight_layout(); save(fig,name)

# c3 detection
fig,ax=plt.subplots(figsize=(3.0,2.6))
groups={"PH":["PH"],"GH":["GH"],"Niche":["TS","MH2","MyoEpi","M"],"Hepato":["F","R","B"],"Other":["PM-like"]}
pct=[(itgb1[a.obs.cell_type.isin(v).values]>0).mean()*100 for v in groups.values()]
ax.bar(range(len(groups)),pct,color="#C44E52",width=0.6)
for i,p_ in enumerate(pct): ax.text(i,p_+0.15,f"{p_:.1f}%",ha="center",fontsize=5.5)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(list(groups),rotation=40); ax.set_ylabel("ITGB1+ cells (%)"); ax.set_title("Receptor detection")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"c"); fig.tight_layout(); save(fig,"c3_ITGB1_detection")

# d dumbbell
fig,ax=plt.subplots(figsize=(4.4,3.0))
piv=gdf.groupby(["source","target","group"])["prob"].sum().unstack(fill_value=0).reset_index()
piv["label"]=piv.source+"→"+piv.target; piv=piv.sort_values("Hemospasia",ascending=True)
for yvar in range(len(piv)): ax.plot([piv.Control.iloc[yvar],piv.Hemospasia.iloc[yvar]],[yvar,yvar],color="#BBB",lw=1.0,zorder=1)
ax.scatter(piv.Control,range(len(piv)),s=28,c="#4C78A8",label="Control",zorder=3)
ax.scatter(piv.Hemospasia,range(len(piv)),s=28,c="#E45756",label="Hemospasia",zorder=3)
ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.label,fontsize=5.8); ax.set_xlabel("Summed ECM–integrin probability")
ax.legend(frameon=False,fontsize=6,loc="lower right"); ax.set_title("Before vs after hemospasia")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"d"); fig.tight_layout(); save(fig,"d_control_vs_hemospasia")

# e1 trajectory tsne
fig,ax=plt.subplots(figsize=(3.7,3.2))
v=t2.obs["projected_pseudotime"].values
sc=ax.scatter(t2.obsm["X_tsne"][:,0],t2.obsm["X_tsne"][:,1],c=v,cmap="viridis",s=1.2,linewidths=0,rasterized=True)
for ct in ["PH","GH","TS","MyoEpi"]:
    m=(t2.obs.cell_type==ct).values; ax.text(np.median(t2.obsm["X_tsne"][m,0]),np.median(t2.obsm["X_tsne"][m,1]),ct,fontsize=6,ha="center",bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
ax.set_xticks([]); ax.set_yticks([]); ax.set_title("Projection onto PH→GH trajectory")
for s in ax.spines.values(): s.set_visible(False)
plt.colorbar(sc,ax=ax,fraction=0.045,pad=0.02).set_label("Projected pseudotime",fontsize=6)
letter(ax,"e"); fig.tight_layout(); save(fig,"e1_trajectory_tsne")

# e2 boxplot
fig,ax=plt.subplots(figsize=(3.2,2.8))
obs=t2.obs
data=[obs.loc[obs.cell_type=="PH","projected_pseudotime"].values,obs.loc[obs.cell_type=="GH","projected_pseudotime"].values,
      obs.loc[obs.compartment=="niche","projected_pseudotime"].values,obs.loc[obs.compartment=="hepatopancreas","projected_pseudotime"].values,
      obs.loc[obs.compartment=="other","projected_pseudotime"].values]
bp=ax.boxplot(data,tick_labels=["PH","GH","Niche","Hepato","Other"],showfliers=False,patch_artist=True,widths=0.6)
for p,c in zip(bp["boxes"],["#E45756","#72B7B2","#4C78A8","#8C510A","#999"]): p.set_facecolor(c); p.set_alpha(0.7)
ax.set_ylabel("Projected pseudotime"); ax.set_title("Compartment comparison"); ax.tick_params(axis="x",rotation=30)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"e"); fig.tight_layout(); save(fig,"e2_compartment_boxplot")
print("done")
