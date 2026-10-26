# -*- coding: utf-8 -*-
"""177_pseudotime_elife_style_try.py
Try an eLife-style PH-rooted pseudotime layout with the final subtype names:
(a) DPT gradient, (b) candidate branch map, (c) DPT distribution.
Exploratory preview only.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import rankdata
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7); sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "elife_try"); os.makedirs(PREV, exist_ok=True)
MAP = {"MH1": "TS", "PH": "PH", "GH": "GH", "MH2": "MH2"}
ORDER = ["PH", "GH", "TS", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}

def normrank(x):
    x=np.asarray(x,float); good=np.isfinite(x); out=np.full(len(x),np.nan)
    r=rankdata(x[good],method="average"); out[good]=(r-1)/max(1,good.sum()-1); return out

def clean(ax,xl,yl):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(xl,fontsize=5.8); ax.set_ylabel(yl,fontsize=5.8)
    for sp in ax.spines.values(): sp.set_color("#777777"); sp.set_linewidth(0.45)

def main():
    A=ad.read_h5ad(os.path.join(OUT,"adata_hemocyte_traj_geodesic.h5ad"))
    A.obs["ct_new"]=pd.Categorical(A.obs["ct"].astype(str).map(MAP),categories=ORDER)
    ct=A.obs["ct_new"].astype(str).values; U=np.asarray(A.obsm["X_umap"])
    root=int(A.uns.get("iroot",0))
    if ct[root]!="PH":
        ph=np.where(ct=="PH")[0]; root=int(ph[np.argmin(((U[ph]-U[ph].mean(0))**2).sum(1))])
    dpt=np.asarray(A.obs["dpt_pseudotime"],float); dpt[~np.isfinite(dpt)]=np.nan
    pt=normrank(dpt); cent={g:U[ct==g].mean(0) for g in ORDER}
    fig,axes=plt.subplots(1,3,figsize=(7.2,3.0),gridspec_kw=dict(wspace=0.32,left=0.06,right=0.98,top=0.88,bottom=0.14))
    # a gradient
    ax=axes[0]; o=np.argsort(-pt)
    s=ax.scatter(U[o,0],U[o,1],c=pt[o],s=1.05,cmap="viridis",linewidths=0,rasterized=True,vmin=0,vmax=1)
    ax.scatter(U[root,0],U[root,1],s=45,marker="*",c="#FF3B30",edgecolors="white",linewidths=0.5,zorder=9)
    ax.annotate("root (PH)",xy=(U[root,0],U[root,1]),xytext=(0,11),textcoords="offset points",fontsize=5.6,color="#B2182B",ha="center",arrowprops=dict(arrowstyle="-",lw=0.5,color="#B2182B"))
    cb=fig.colorbar(s,ax=ax,fraction=0.046,pad=0.03); cb.set_label("DPT (rank display)",fontsize=5.2); cb.ax.tick_params(labelsize=4.8,length=0.9); cb.outline.set_visible(False)
    clean(ax,"UMAP 1","UMAP 2"); ax.set_title("PH-rooted pseudotime",fontsize=6.8,pad=4)
    # b group map + candidate paths
    ax=axes[1]
    for g in ORDER:
        m=ct==g; ax.scatter(U[m,0],U[m,1],s=1.05,c=COL[g],linewidths=0,rasterized=True,label="%s (%d)"%(g,m.sum()))
    paths=[(cent["PH"],cent["TS"],-0.12,"#2B2B2B",1.1),(cent["TS"],cent["GH"],0.30,"#2166AC",0.95),(cent["TS"],cent["MH2"],-0.30,"#8C510A",0.95)]
    for p1,p2,rad,col,lw in paths:
        ax.add_patch(FancyArrowPatch(p1,p2,connectionstyle="arc3,rad=%s"%rad,arrowstyle="-|>,head_length=0.45,head_width=0.20",mutation_scale=8,lw=lw,color=col,alpha=0.9,zorder=10,shrinkA=2,shrinkB=2,path_effects=[pe.withStroke(linewidth=2.0,foreground="white")]))
    offs={"PH":(-12,9),"GH":(0,0),"TS":(0,0),"MH2":(0,0)}
    for g in ORDER:
        c=cent[g]; ax.annotate(g,xy=(c[0],c[1]),xytext=offs[g],textcoords="offset points",fontsize=5.9,fontweight="bold",ha="center",va="center",color="white",path_effects=[pe.withStroke(linewidth=1.9,foreground="#333333")],zorder=11)
    ax.legend(loc="lower left",bbox_to_anchor=(-0.01,-0.02),ncol=2,frameon=False,markerscale=4.0,fontsize=4.9,handletextpad=0.2,columnspacing=0.6,labelspacing=0.25)
    clean(ax,"UMAP 1",""); ax.set_title("Candidate DPT paths",fontsize=6.8,pad=4)
    # c distribution
    ax=axes[2]
    data=[pt[ct==g] for g in ORDER]
    vp=ax.violinplot(data,positions=np.arange(1,5),widths=0.68,showextrema=False,showmeans=False)
    for body,g in zip(vp["bodies"],ORDER): body.set_facecolor(COL[g]); body.set_edgecolor("none"); body.set_alpha(0.55)
    rng=np.random.default_rng(0)
    for i,g in enumerate(ORDER,start=1):
        d=data[i-1]; j=rng.uniform(-0.10,0.10,size=len(d))
        ax.scatter(i+j,d,s=0.8,c=COL[g],alpha=0.35,linewidths=0,rasterized=True)
        ax.scatter([i],[np.nanmedian(d)],s=15,c="white",edgecolors=COL[g],linewidths=0.6,zorder=5)
    ax.set_xticks(range(1,5)); ax.set_xticklabels(ORDER,fontsize=5.2)
    ax.set_ylim(0,1); ax.set_ylabel("DPT (rank display)",fontsize=5.6); ax.set_xlabel("")
    ax.tick_params(axis="y",labelsize=4.9,length=1.1); ax.tick_params(axis="x",length=0)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.set_title("Pseudotime by state",fontsize=6.8,pad=4)
    for ax,lab in zip(axes,["a","b","c"]):
        ax.text(-0.18,1.03,lab,transform=ax.transAxes,fontsize=8,fontweight="bold",va="top")
    fig.text(0.5,0.018,"Candidate paths are based on DPT ordering; RNA velocity supports PH as source, not the downstream branch order.",ha="center",va="bottom",fontsize=5.0,color="#333333")
    p=os.path.join(PREV,"Fig_pseudotime_elife_style_try.png"); fig.savefig(p,dpi=350,bbox_inches=None,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__":
    log=open(os.path.join(r"D:\桌面\single cell\logs","pseudotime_elife_style_try.log"),"w",encoding="utf8",buffering=1)
    sys.stdout=log; sys.stderr=log; main(); log.close()
