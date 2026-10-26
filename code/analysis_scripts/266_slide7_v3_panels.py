# -*- coding: utf-8 -*-
"""266_slide7_v3_panels.py — standalone v3 panels a-g."""
import os
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1")
BASE=os.path.join(OUT,"figures_pub","Slide7_v3","panels")
for s in ["svg","tiff","png"]: os.makedirs(os.path.join(BASE,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def save(fig,n):
    fig.savefig(os.path.join(BASE,"svg",n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(BASE,"tiff",n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(BASE,"png",n+".png"),dpi=220,bbox_inches="tight"); plt.close(fig); print("saved",n,flush=True)
def letter(ax,s): ax.text(-0.08,1.05,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

t=ad.read_h5ad(os.path.join(OUT,"adata_hemocyte_traj_geodesic.h5ad"))
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv")); ph.columns=["cell","group","subcluster","umap1","umap2"]
meta=ph.set_index("cell")[["group","subcluster"]].rename(columns={"group":"group_ph"}); t.obs=t.obs.join(meta)
phd=t.obs[t.obs.ct=="PH"].copy(); phd["group"]=phd["group_ph"]; phd["subcluster"]=pd.to_numeric(phd["subcluster"])

# a
fig,axes=plt.subplots(2,2,figsize=(3.6,3.1))
for k,ax in enumerate(axes.ravel()):
    for grp,col in [("Control","#4C78A8"),("Hemospasia","#E45756")]:
        v=phd.loc[(phd.subcluster==k)&(phd.group==grp),"dpt_pseudotime"].dropna().values
        if len(v)>10:
            x=np.linspace(max(0,v.min()-0.01),min(1,v.max()+0.01),120)
            ax.plot(x,gaussian_kde(v)(x),color=col,lw=1.0); ax.axvline(np.median(v),color=col,lw=0.5,ls=":")
    ax.set_title(f"PH_s{k}",fontsize=6.5); ax.set_yticks([]); ax.set_xlabel("DPT",labelpad=1)
    if k%2==0: ax.set_ylabel("Density",labelpad=1)
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    if k==0: ax.legend(handles=[Patch(color="#4C78A8",label="Control"),Patch(color="#E45756",label="Hemospasia")],frameon=False,fontsize=5,loc="upper right")
letter(axes[0,0],"a"); fig.tight_layout(); save(fig,"panel_a_PH_density")

# b
h=pd.read_csv(os.path.join(OUT,"Fig5_PH_dpt_heatmap_data.csv")); piv=h.pivot_table(index="gene",columns="bin",values="zscore")
piv=piv.reindex(h.groupby("gene")["median_dpt"].mean().sort_values().index.tolist())
fig,ax=plt.subplots(figsize=(3.4,3.2)); im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index,fontsize=5.2)
ax.set_xticks([0,len(piv.columns)-1]); ax.set_xticklabels(["early","late"]); ax.set_xlabel("Pseudotime")
ax.annotate("",xy=(len(piv.columns)-0.5,-1.0),xytext=(-0.5,-1.0),arrowprops=dict(arrowstyle="->",lw=0.9),annotation_clip=False)
plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02).set_label("z",fontsize=5.5); letter(ax,"b"); fig.tight_layout(); save(fig,"panel_b_PH_heatmap")

# c grouped
cc=pd.read_csv(os.path.join(OUT,"Fig5_dynamic_cellcycle_by_group.csv"))
fig,ax=plt.subplots(figsize=(4.0,2.6)); pops=["PH","MH2","R","PM-like"]; xpos=[]; labs=[]; pos=0
for p_ in pops:
    for g_ in ["Control","Hemospasia"]:
        r=cc[(cc.population==p_)&(cc.group==g_)].iloc[0]; left=0
        for phas,col in zip(["pct_G1","pct_S","pct_G2M"],["#BAB0AC","#72B7B2","#E45756"]):
            ax.bar(pos,r[phas],bottom=left,color=col,width=0.72,edgecolor="white",linewidth=0.3); left+=r[phas]
        xpos.append(pos); labs.append(f"{p_}\n{'C' if g_=='Control' else 'H'}"); pos+=1
    pos+=0.35
ax.set_xticks(xpos); ax.set_xticklabels(labs,fontsize=5.0); ax.set_ylim(0,100); ax.set_ylabel("Cells (%)")
ax.legend(handles=[Patch(color="#BAB0AC",label="G1"),Patch(color="#72B7B2",label="S"),Patch(color="#E45756",label="G2M")],frameon=False,fontsize=5,ncol=3,loc="upper center",bbox_to_anchor=(0.5,1.20))
for i in range(0,len(xpos),2):
    d=cc[(cc.population==pops[i//2])&(cc.group=="Hemospasia")].iloc[0].pct_G2M-cc[(cc.population==pops[i//2])&(cc.group=="Control")].iloc[0].pct_G2M
    ax.text((xpos[i]+xpos[i+1])/2,102,f"Δ{d:+.1f}",ha="center",fontsize=5.2,color="#C44E52")
ax.text(0.5,-0.22,"C = Control; H = Hemospasia",transform=ax.transAxes,ha="center",fontsize=4.8)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.text(-0.16,1.10,"c",transform=ax.transAxes,fontsize=10,fontweight="bold",va="top"); fig.tight_layout(); save(fig,"panel_c_cellcycle_group")

# d
me=pd.read_csv(os.path.join(OUT,"Fig5_coexpression_module_effects.csv")); piv=me.pivot_table(index="parent",columns="module",values="cohens_d")
mod={"M1":"M1 PM-like/cuticle","M2":"M2 ECM/lipid","M3":"M3 Absorptive","M4":"M4 Mixed","M5":"M5 ECM/contractile","M6":"M6 Transport","M7":"M7 PH progenitor"}
fig,ax=plt.subplots(figsize=(4.4,2.6)); im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-1.5,vmax=1.5)
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels([mod.get(c,c) for c in piv.columns],rotation=25,ha="right",fontsize=5.2)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v=piv.values[i,j]
        if np.isfinite(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=4.8,color="white" if abs(v)>0.9 else "#333")
plt.colorbar(im,ax=ax,fraction=0.03,pad=0.015).set_label("Cohen's d",fontsize=5.5); letter(ax,"d"); fig.tight_layout(); save(fig,"panel_d_modules_named")

# e UpSet
deg=pd.read_csv(os.path.join(OUT,"Fig5_DEG_shared_unique.csv")); cnt=deg.groupby("combination").size().sort_values(ascending=False).head(8)
sets=["PH","MH2","R","PM-like"]; combos=cnt.index.tolist(); vals=cnt.values
fig=plt.figure(figsize=(3.8,3.0)); gs=fig.add_gridspec(2,1,height_ratios=[1,1.55],hspace=0.06)
ax1=fig.add_subplot(gs[0]); ax1.bar(range(len(combos)),vals,color="#4C78A8",width=0.62)
for i,v in enumerate(vals): ax1.text(i,v+1,str(v),ha="center",fontsize=5)
ax1.set_xticks([]); ax1.set_ylabel("DEGs"); ax1.set_ylim(0,max(vals)*1.25)
for s in ["top","right"]: ax1.spines[s].set_visible(False)
ax2=fig.add_subplot(gs[1])
for j,c in enumerate(combos):
    mem=[s for s in sets if s in c.split("+")]; ys=[sets.index(s) for s in mem]
    ax2.scatter([j]*len(ys),ys,s=18,c="#C44E52",zorder=3)
    if len(ys)>1: ax2.plot([j,j],[min(ys),max(ys)],color="#C44E52",lw=0.7,zorder=2)
    for s in sets:
        if s not in mem: ax2.scatter([j],[sets.index(s)],s=10,c="#DDDDDD",zorder=3)
short={"PH+MH2+R+PM-like":"All 4","PM-like":"PM-like","PH+MH2+R":"PH+MH2+R","PH":"PH","R":"R","PH+MH2":"PH+MH2","MH2":"MH2","PH+MH2+PM-like":"PH+MH2+PM"}
ax2.set_yticks(range(len(sets))); ax2.set_yticklabels(sets); ax2.set_xticks(range(len(combos)))
ax2.set_xticklabels([short.get(c,c) for c in combos],rotation=35,ha="right",fontsize=4.8)
ax2.set_xlim(-0.6,len(combos)-0.4); ax2.invert_yaxis()
for s in ["top","right"]: ax2.spines[s].set_visible(False)
letter(ax1,"e"); fig.tight_layout(); save(fig,"panel_e_DEG_upset")

# f
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv")); agg=lr.groupby(["source","target"],as_index=False)["prob"].sum()
srcs=["TS","MH2","MyoEpi","M"]; tgts=["PH","GH"]; cols={"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2"}
s_y={x:y for x,y in zip(srcs,np.linspace(1,0,len(srcs)))}; t_y={x:y for x,y in zip(tgts,np.linspace(0.75,0.25,len(tgts)))}
fig,ax=plt.subplots(figsize=(3.2,3.0)); ax.set_xlim(0,1); ax.set_ylim(-0.1,1.1); ax.axis("off")
for s in srcs: ax.scatter(0.12,s_y[s],s=130,color=cols[s],zorder=3); ax.text(0.08,s_y[s],s,ha="right",va="center",fontsize=6.5)
for t_ in tgts: ax.scatter(0.88,t_y[t_],s=150,color="#E45756" if t_=="PH" else "#72B7B2",zorder=3); ax.text(0.92,t_y[t_],t_,ha="left",va="center",fontsize=6.5)
for r in agg.itertuples():
    if r.prob<0.01: continue
    p0=np.array([0.14,s_y[r.source]]); p1=np.array([0.86,t_y[r.target]]); ctrl=np.array([0.5,(p0[1]+p1[1])/2+0.06])
    ax.add_patch(PathPatch(Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3]),fill=False,color=cols[r.source],lw=0.4+r.prob*7,alpha=0.65))
ax.text(0.5,1.06,"Niche cells",ha="center",fontsize=7,fontweight="bold"); ax.text(0.5,0.98,"→ Hemocytes",ha="center",fontsize=6.5); letter(ax,"f"); fig.tight_layout(); save(fig,"panel_f_LR_bipartite")

# g alluvial with p-value table (no overlapping inline labels)
PROGS=["Progenitor","Macrophage-like","Hemocytes","ECM-producing","Capsular","Epithelial","SLC-rich","Lymphocyte precursor-like"]
ORGS=["Hepatic sinus","Hemolymph","Hepatopancreas"]; O2=["Hs","Hemolymph","Hp"]
COL={"Hepatic sinus":"#B2182B","Hemolymph":"#2166AC","Hepatopancreas":"#8C510A"}
CELL={"Hepatic sinus":{"PH":"PH","MH2":"MH2","R":"R","PM-like":"PM-like","MyoEpi":"MyoEpi"},"Hemolymph":{"1":"GH-like","6":"MH1-like","11":"cycling","10":"c10","9":"c9"},"Hepatopancreas":{"6":"prolif.","2":"epi/immune","7":"R-like","15":"neuroendo.","14":"c14"}}
R=pd.read_csv(os.path.join(OUT,"Fig6_program_organ_summary.csv")); P=pd.read_csv(os.path.join(OUT,"Fig6_permutation_results.csv")); S=pd.read_csv(os.path.join(OUT,"Fig6_marker_sensitivity.csv"))
def bez(x0,y0,x1,y1):
    t=np.linspace(0,1,80); x=(1-t)**3*x0+3*(1-t)**2*t*(x0+0.28*(x1-x0))+3*(1-t)*t**2*(x1-0.28*(x1-x0))+t**3*x1
    y=(1-t)**3*y0+3*(1-t)**2*t*y0+3*(1-t)*t**2*y1+t**3*y1; return x,y
def stab(p,o):
    d=S[(S.program==p)&(S.organ==o)].dropna(subset=["best_group"]); b=d.best_group.astype(str).tolist(); return len(b)>0 and len(set(b))==1
fig,ax=plt.subplots(figsize=(5.0,3.3)); ax.axis("off"); ax.set_xlim(0,1); ax.set_ylim(0,1)
ly={p:0.90-0.105*i for i,p in enumerate(PROGS)}; my={"Hepatic sinus":0.78,"Hemolymph":0.50,"Hepatopancreas":0.22}
for p,y in ly.items(): ax.scatter([0.06],[y],s=18,c="#555",zorder=8); ax.text(0.048,y,p,ha="right",va="center",fontsize=4.4)
for o,y in my.items(): ax.scatter([0.26],[y],s=32,c=COL[o],zorder=8); ax.text(0.26,y+0.055,O2[ORGS.index(o)],ha="center",va="bottom",fontsize=4.4,color=COL[o],weight="bold")
for p in PROGS:
    a=P[(P.program==p)&(P.organ.isin(ORGS))].dropna(subset=["score"]).sort_values("score",ascending=False)
    sig=a[a.q_perm<0.10]; sub=sig.copy() if len(sig)>=2 else pd.concat([sig,a.head(2)],ignore_index=True).drop_duplicates(subset=["organ"]).sort_values("score",ascending=False)
    for r in sub.itertuples():
        x,y=bez(0.06,ly[p],0.26,my[r.organ]); st=stab(p,r.organ); q=float(r.q_perm) if pd.notna(r.q_perm) else 1.0
        if not st: ls=(0,(1,2)); lw=0.7; alpha=0.18
        elif q<0.05: ls="-"; lw=1.6+1.1*max(0,float(r.score))/4; alpha=0.90
        elif q<0.10: ls=(0,(4,2)); lw=1.2; alpha=0.75
        else: ls="-"; lw=0.7; alpha=0.25
        ax.plot(x,y,lw=lw,ls=ls,color=COL[r.organ],alpha=alpha,zorder=2,solid_capstyle="round")
for o in ORGS:
    sub=R[R.organ==o].dropna(subset=["z_score"]).copy(); sub["label"]=sub.best_cluster.astype(str).map(CELL[o]).fillna(sub.best_cluster.astype(str))
    ag=sub.groupby("label",as_index=False).agg(score=("z_score","sum")).sort_values("score",ascending=False)
    ys=np.linspace(my[o]+0.07,my[o]-0.07,len(ag)) if len(ag)>1 else [my[o]]
    for (lab,sc),yy in zip(ag.itertuples(index=False),ys):
        x,y=bez(0.26,my[o],0.44,yy); ax.plot(x,y,lw=0.5+2.4*float(sc)/8,color=COL[o],alpha=0.5,zorder=3)
        ax.scatter([0.44],[yy],s=10,c=COL[o],edgecolors="white",linewidths=0.25,zorder=8); ax.text(0.452,yy,lab,ha="left",va="center",fontsize=3.7,color=COL[o])
# p-value table
ax.text(0.56,0.965,"Permutation P",fontsize=5.2,fontweight="bold")
for j,o in enumerate(ORGS): ax.text(0.58+0.075*j,0.915,O2[j],fontsize=4.2,color=COL[o],ha="center")
for p in PROGS:
    y=ly[p]
    for j,o in enumerate(ORGS):
        r=P[(P.program==p)&(P.organ==o)]
        txt="—"
        if len(r) and pd.notna(r.iloc[0].q_perm) and r.iloc[0].q_perm<0.10:
            txt=f"{r.iloc[0].p_perm:.3g}"
        ax.text(0.58+0.075*j,y,txt,fontsize=3.5,ha="center",va="center",color=COL[o] if txt!="—" else "#AAAAAA")
ax.legend(handles=[Line2D([0],[0],color="#333",lw=1.6,label="FDR<0.05"),Line2D([0],[0],color="#333",lw=1.2,ls=(0,(4,2)),label="trend 0.05-0.10")],frameon=False,fontsize=3.6,loc="lower left",bbox_to_anchor=(0.0,-0.02),handlelength=1.2)
letter(ax,"g"); fig.tight_layout(); save(fig,"panel_g_LyO_alluvial")
print("done ->",BASE)
