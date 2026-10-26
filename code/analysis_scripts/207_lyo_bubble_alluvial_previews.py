# -*- coding: utf-8 -*-
"""207_lyo_bubble_alluvial_previews.py
Two intuitive views of the LyO-program cross-organ summary.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
R=pd.read_csv(os.path.join(OUT,"Fig6_program_organ_summary.csv"))
PROGS=["Progenitor","Macrophage-like","Hemocytes","ECM-producing","Capsular","Epithelial","SLC-rich","Lymphocyte precursor-like"]
ORGS=["Hepatic sinus","Hemolymph","Hepatopancreas"]
COL={"Hepatic sinus":"#B2182B","Hemolymph":"#2166AC","Hepatopancreas":"#8C510A"}

def bubble():
 fig=plt.figure(figsize=(4.8,3.4),facecolor="white"); ax=fig.add_axes([0.24,0.17,0.66,0.72])
 for i,p in enumerate(PROGS):
  for j,o in enumerate(ORGS):
   r=R[(R.program==p)&(R.organ==o)]
   if len(r)==0 or pd.isna(r.iloc[0].z_score):
    ax.scatter(j,i,s=12,facecolors="none",edgecolors="#BBBBBB",linewidths=0.6)
    ax.text(j,i,"n.d.",ha="center",va="center",fontsize=4.1,color="#999999")
   else:
    r=r.iloc[0]; z=float(r.z_score); sz=8+12*z
    ax.scatter(j,i,s=sz,c=z,cmap="viridis",vmin=0,vmax=4,edgecolors="white",linewidths=0.5,zorder=3)
    lab=str(r.best_cluster)
    ax.text(j,i,lab,ha="center",va="center",fontsize=4.4,color="white" if z<2.4 else "#222222",zorder=4)
 ax.set_xticks(range(len(ORGS))); ax.set_xticklabels(["Hepatic\nsinus","Hemolymph","Hepato-\npancreas"],fontsize=5.2)
 ax.set_yticks(range(len(PROGS))); ax.set_yticklabels(PROGS,fontsize=5.1)
 ax.set_xlim(-0.55,2.55); ax.set_ylim(len(PROGS)-0.5,-0.5)
 ax.tick_params(length=0)
 for sp in ax.spines.values(): sp.set_visible(False)
 ax.grid(color="#EEEEEE",lw=0.4); ax.set_axisbelow(True)
 ax.set_title("LyO programmes across three organs",fontsize=6.5,pad=6)
 cb=fig.colorbar(plt.cm.ScalarMappable(cmap="viridis",norm=plt.Normalize(0,4)),ax=ax,fraction=0.035,pad=0.03)
 cb.set_label("best-match z-score",fontsize=4.6); cb.ax.tick_params(labelsize=4.2,length=0.7); cb.outline.set_visible(False)
 fig.text(0.50,0.025,"Bubble label: best matching cluster. LyO markers are the reference; n.d. = no matched genes.",ha="center",fontsize=4.5,color="#555555")
 p=os.path.join(PREV,"Fig6_bubble_matrix_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

def alluvial():
 fig=plt.figure(figsize=(5.2,4.0),facecolor="white"); ax=fig.add_axes([0.02,0.08,0.96,0.82]); ax.axis("off")
 ly={p:1-0.11*(i+0.5) for i,p in enumerate(PROGS)}
 oy={"Hepatic sinus":0.78,"Hemolymph":0.50,"Hepatopancreas":0.22}
 # background nodes
 for p,y in ly.items():
  ax.text(0.055,y,p,ha="right",va="center",fontsize=4.8)
  ax.scatter([0.075],[y],s=22,c="#555555",zorder=5)
 for o,y in oy.items():
  ax.text(0.925,y,o,ha="left",va="center",fontsize=5.0,weight="bold",color=COL[o])
  ax.scatter([0.875],[y],s=30,c=COL[o],zorder=5)
 # top-2 links per program
 for p in PROGS:
  sub=R[R.program==p].dropna(subset=["z_score"]).sort_values("z_score",ascending=False).head(2)
  for r in sub.itertuples():
   y0=ly[p]; y1=oy[r.organ]
   # cubic bezier
   verts=[(0.075,y0),(0.38,y0),(0.60,y1),(0.875,y1)]
   codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
   patch=PathPatch(Path(verts,codes),fill=False,lw=0.6+3.8*max(0,float(r.z_score))/4.0,
                   color=COL[r.organ],alpha=0.45,zorder=2)
   ax.add_patch(patch)
 ax.set_xlim(0,1); ax.set_ylim(0,1)
 ax.set_title("Shared programme flows: LyO → target organs",fontsize=6.5,pad=4)
 fig.legend(handles=[Patch(color=COL[o],label=o) for o in ORGS],loc="lower center",bbox_to_anchor=(0.5,0.015),ncol=3,frameon=False,fontsize=4.8,handlelength=0.8,columnspacing=0.8)
 fig.text(0.5,0.045,"Top 2 organ matches per LyO programme; link width = z-score. Marker-level candidate sharing, not migration.",ha="center",fontsize=4.4,color="#555555")
 p=os.path.join(PREV,"Fig6_alluvial_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": bubble(); alluvial()
