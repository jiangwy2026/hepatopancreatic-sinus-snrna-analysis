# -*- coding: utf-8 -*-
"""208_three_column_alluvial.py
Three-column alluvial: LyO programme -> organ -> matched cell state.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from PIL import Image
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
R=pd.read_csv(os.path.join(OUT,"Fig6_program_organ_summary.csv"))
PROGS=["Progenitor","Macrophage-like","Hemocytes","ECM-producing","Capsular","Epithelial","SLC-rich","Lymphocyte precursor-like"]
ORGS=["Hepatic sinus","Hemolymph","Hepatopancreas"]
O2=["Hepatopancreatic sinus (Hs)","Hemolymph","Hepatopancreas"]
COL={"Hepatic sinus":"#B2182B","Hemolymph":"#2166AC","Hepatopancreas":"#8C510A"}
CELL={
 "Hepatic sinus":{"PH":"PH","MH2":"MH2","R":"R","PM-like":"PM-like","MyoEpi":"MyoEpi"},
 "Hemolymph":{"1":"GH-like (c1)","6":"MH1-like (c6)","11":"cycling/PH-like (c11)","10":"c10 (unknown)","9":"c9 (unknown)"},
 "Hepatopancreas":{"6":"proliferative (c6)","2":"epithelial/immune (c2)","7":"R-like (c7)","15":"neuroendocrine (c15)","14":"c14 (unknown)"},
}

def bezier(x0,y0,x1,y1):
 t=np.linspace(0,1,80); x=(1-t)**3*x0+3*(1-t)**2*t*(x0+0.35*(x1-x0))+3*(1-t)*t**2*(x1-0.35*(x1-x0))+t**3*x1
 y=(1-t)**3*y0+3*(1-t)**2*t*y0+3*(1-t)*t**2*y1+t**3*y1
 return x,y

def main():
 fig=plt.figure(figsize=(5.4,4.05),facecolor="white"); ax=fig.add_axes([0.02,0.10,0.96,0.80]); ax.axis("off")
 ax.set_xlim(0,1); ax.set_ylim(0,1)
 left_y={p:0.88-0.108*i for i,p in enumerate(PROGS)}
 mid_y={"Hepatic sinus":0.78,"Hemolymph":0.50,"Hepatopancreas":0.22}
 # left nodes and labels
 for p,y in left_y.items():
  ax.scatter([0.075],[y],s=26,c="#555555",zorder=8)
  ax.text(0.06,y,p,ha="right",va="center",fontsize=4.8)
 # middle nodes and organ labels
 for o,y in mid_y.items():
  col=COL[o]
  ax.scatter([0.30],[y],s=45,c=col,zorder=8)
  ax.text(0.30,y+0.055,O2[ORGS.index(o)],ha="center",va="bottom",fontsize=4.7,color=col,weight="bold")
 # top-2 program->organ links
 for p in PROGS:
  sub=R[R.program==p].dropna(subset=["z_score"]).sort_values("z_score",ascending=False).head(2)
  for r in sub.itertuples():
   x,y=bezier(0.075,left_y[p],0.30,mid_y[r.organ])
   ax.plot(x,y,lw=0.5+3.4*max(0,float(r.z_score))/4.0,color=COL[r.organ],alpha=0.40,zorder=2,solid_capstyle="round")
 # right nodes: aggregate by organ + best_cluster
 for o in ORGS:
  sub=R[(R.organ==o)].dropna(subset=["z_score"]).copy()
  sub["label"]=sub.best_cluster.astype(str).map(CELL[o]).fillna(sub.best_cluster.astype(str)+" (unknown)")
  agg=sub.groupby("label",as_index=False).agg(score=("z_score","sum"))
  agg=agg.sort_values("score",ascending=False)
  ys=np.linspace(mid_y[o]+0.075,mid_y[o]-0.075,len(agg)) if len(agg)>1 else [mid_y[o]]
  for (lab,sc),yy in zip(agg.itertuples(index=False),ys):
   y0=mid_y[o]; x,y=bezier(0.30,y0,0.51,yy)
   ax.plot(x,y,lw=0.6+4.0*float(sc)/8.0,color=COL[o],alpha=0.55,zorder=3,solid_capstyle="round")
   ax.scatter([0.51],[yy],s=18,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8)
   ax.text(0.525,yy,lab,ha="left",va="center",fontsize=4.5,color=COL[o])
 ax.set_title("Shared programme flows: LyO → organ → matched cell state",fontsize=6.6,pad=4)
 fig.legend(handles=[Patch(color=COL[o],label=O2[i]) for i,o in enumerate(ORGS)],loc="lower center",bbox_to_anchor=(0.5,0.010),ncol=3,frameon=False,fontsize=4.8,handlelength=0.8,columnspacing=0.8)
 fig.text(0.5,0.075,"Top 2 organ matches per LyO programme. Link width = z-score; marker-level candidate sharing, not migration.",ha="center",fontsize=4.4,color="#555555")
 p=os.path.join(PREV,"Fig6_alluvial_threecolumn_compact_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); print("saved",p)
 deliv=os.path.join(r"D:\桌面\single cell\文章配图", "最终_2026-10-19", "四器官模型_河流图"); os.makedirs(deliv,exist_ok=True)
 stem=os.path.join(deliv,"Fig6_alluvial_threecolumn_compact")
 fig.savefig(stem+".svg",bbox_inches="tight",pad_inches=0.04,facecolor="white")
 fig.savefig(stem+".tiff",dpi=1200,bbox_inches="tight",pad_inches=0.04,facecolor="white",pil_kwargs={"compression":"tiff_lzw"})
 with Image.open(stem+".tiff") as im:
  if "A" in im.getbands():
   rgba=im.convert("RGBA"); bg=Image.new("RGB",rgba.size,(255,255,255)); bg.paste(rgba,mask=rgba.getchannel("A")); tmp=stem+".tmp.tiff"; bg.save(tmp,format="TIFF",dpi=(1200,1200),compression="tiff_lzw"); os.replace(tmp,stem+".tiff")
 print("saved",stem+".svg"); print("saved",stem+".tiff")
 plt.close(fig)

if __name__=="__main__": main()
