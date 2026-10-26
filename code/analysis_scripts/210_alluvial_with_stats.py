# -*- coding: utf-8 -*-
"""210_alluvial_with_stats.py
Three-column alluvial with permutation p-values, significance and stability encoding.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
R=pd.read_csv(os.path.join(OUT,"Fig6_program_organ_summary.csv"))
P=pd.read_csv(os.path.join(OUT,"Fig6_permutation_results.csv"))
S=pd.read_csv(os.path.join(OUT,"Fig6_marker_sensitivity.csv"))
PROGS=["Progenitor","Macrophage-like","Hemocytes","ECM-producing","Capsular","Epithelial","SLC-rich","Lymphocyte precursor-like"]
ORGS=["Hepatic sinus","Hemolymph","Hepatopancreas"]
O2=["Hepatopancreatic sinus (Hs)","Hemolymph","Hepatopancreas"]
COL={"Hepatic sinus":"#B2182B","Hemolymph":"#2166AC","Hepatopancreas":"#8C510A"}
CELL={"Hepatic sinus":{"PH":"PH","MH2":"MH2","R":"R","PM-like":"PM-like","MyoEpi":"MyoEpi"},"Hemolymph":{"1":"GH-like (c1)","6":"MH1-like (c6)","11":"cycling/PH-like (c11)","10":"c10 (unknown)","9":"c9 (unknown)"},"Hepatopancreas":{"6":"proliferative (c6)","2":"epithelial/immune (c2)","7":"R-like (c7)","15":"neuroendocrine (c15)","14":"c14 (unknown)"}}

def bezier(x0,y0,x1,y1):
 t=np.linspace(0,1,80); x=(1-t)**3*x0+3*(1-t)**2*t*(x0+0.28*(x1-x0))+3*(1-t)*t**2*(x1-0.28*(x1-x0))+t**3*x1
 y=(1-t)**3*y0+3*(1-t)**2*t*y0+3*(1-t)*t**2*y1+t**3*y1
 return x,y

def stability(p,o):
 d=S[(S.program==p)&(S.organ==o)].dropna(subset=["best_group"])
 b=d.best_group.astype(str).tolist()
 return (len(b)>0 and len(set(b))==1)

def main():
 fig=plt.figure(figsize=(5.8,4.25),facecolor="white"); ax=fig.add_axes([0.02,0.11,0.96,0.79]); ax.axis("off"); ax.set_xlim(0,1); ax.set_ylim(0,1)
 left_y={p:0.88-0.108*i for i,p in enumerate(PROGS)}; mid_y={"Hepatic sinus":0.78,"Hemolymph":0.50,"Hepatopancreas":0.22}
 for p,y in left_y.items():
  ax.scatter([0.075],[y],s=24,c="#555555",zorder=8); ax.text(0.06,y,p,ha="right",va="center",fontsize=4.7)
 for o,y in mid_y.items():
  ax.scatter([0.30],[y],s=42,c=COL[o],zorder=8); ax.text(0.30,y+0.055,O2[ORGS.index(o)],ha="center",va="bottom",fontsize=4.6,color=COL[o],weight="bold")
 # left->middle links with stats
 for p in PROGS:
  allsub=P[(P.program==p)&(P.organ.isin(ORGS))].dropna(subset=["score"]).sort_values("score",ascending=False)
  sig=allsub[allsub.q_perm<0.10]
  if len(sig)>=2:
   sub=sig.copy()
  else:
   sub=pd.concat([sig,allsub.head(2)],ignore_index=True).drop_duplicates(subset=["organ"]).sort_values("score",ascending=False)
  nlinks=len(sub)
  for k,r in enumerate(sub.itertuples()):
   x,y=bezier(0.075,left_y[p],0.30,mid_y[r.organ])
   st=stability(p,r.organ); q=float(r.q_perm) if pd.notna(r.q_perm) else 1.0
   if not st:
    ls=(0,(1,2)); lw=0.8; alpha=0.18
   elif q<0.05:
    ls="-"; lw=2.1+1.5*max(0,float(r.score))/4; alpha=0.92
   elif q<0.10:
    ls=(0,(4,2)); lw=1.5; alpha=0.75
   else:
    ls="-"; lw=0.8; alpha=0.25
   ax.plot(x,y,lw=lw,ls=ls,color=COL[r.organ],alpha=alpha,zorder=2,solid_capstyle="round")
   # p label near the left node, stacked per programme to avoid mid-link collisions
   yy=left_y[p]+(k-(nlinks-1)/2)*0.018
   ax.text(0.115,yy,"p=%.3g"%r.p_perm,ha="left",va="center",fontsize=4.0,color=COL[r.organ],
           bbox=dict(boxstyle="round,pad=0.10",fc="white",ec="none",alpha=0.72),zorder=9)
 # middle->right aggregate
 for o in ORGS:
  sub=R[(R.organ==o)].dropna(subset=["z_score"]).copy(); sub["label"]=sub.best_cluster.astype(str).map(CELL[o]).fillna(sub.best_cluster.astype(str)+" (unknown)")
  agg=sub.groupby("label",as_index=False).agg(score=("z_score","sum")).sort_values("score",ascending=False)
  ys=np.linspace(mid_y[o]+0.075,mid_y[o]-0.075,len(agg)) if len(agg)>1 else [mid_y[o]]
  for (lab,sc),yy in zip(agg.itertuples(index=False),ys):
   x,y=bezier(0.30,mid_y[o],0.51,yy)
   ax.plot(x,y,lw=0.6+3.2*float(sc)/8,color=COL[o],alpha=0.50,zorder=3,solid_capstyle="round")
   ax.scatter([0.51],[yy],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8)
   ax.text(0.525,yy,lab,ha="left",va="center",fontsize=4.4,color=COL[o])
 ax.set_title("LyO programme sharing with statistical support",fontsize=6.5,pad=4)
 # legends
 organ_handles=[Patch(color=COL[o],label=O2[i]) for i,o in enumerate(ORGS)]
 style_handles=[Line2D([0],[0],color="#333333",lw=2.2,label="FDR < 0.05"),
                Line2D([0],[0],color="#333333",lw=1.5,ls=(0,(4,2)),label="trend 0.05–0.10"),
                Line2D([0],[0],color="#999999",lw=0.9,ls=(0,(1,2)),alpha=0.35,label="unstable"),
                Line2D([0],[0],color="#777777",lw=0.8,alpha=0.25,label="stable n.s.")]
 fig.legend(handles=organ_handles,loc="lower center",bbox_to_anchor=(0.5,0.012),ncol=3,frameon=False,fontsize=4.6,handlelength=0.8,columnspacing=0.7,title="Organ",title_fontsize=4.4)
 fig.legend(handles=style_handles,loc="lower left",bbox_to_anchor=(0.03,0.012),ncol=1,frameon=False,fontsize=4.2,handlelength=1.5,labelspacing=0.25)
 fig.text(0.50,0.075,"Numbers: permutation p. All significant/trend links plus top-2 organ matches are shown; marker-level candidate sharing, not migration.",ha="center",fontsize=4.2,color="#555555")
 p=os.path.join(PREV,"Fig6_alluvial_threecolumn_stats_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": main()
