# -*- coding: utf-8 -*-
"""211_alluvial_supported_only.py
Only statistically supported/trend LyO programme-organ links are drawn.
Both left and right segments of each flow share the same style.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D
from PIL import Image
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
mpl.rcParams.update({
 "font.family": "serif",
 "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
 "svg.fonttype": "none",
 "pdf.fonttype": 42,
})
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
P=pd.read_csv(os.path.join(OUT,"Fig6_permutation_results.csv"))
R=pd.read_csv(os.path.join(OUT,"Fig6_program_organ_summary.csv"))
ORGS=["Hepatic sinus","Hemolymph","Hepatopancreas"]
O2=["Hepatopancreatic sinus (Hs)","Hemolymph","Hepatopancreas"]
COL={"Hepatic sinus":"#B2182B","Hemolymph":"#2166AC","Hepatopancreas":"#8C510A"}
CELL={"Hepatic sinus":{"PH":"PH","MH2":"MH2","R":"R","PM-like":"PM-like","MyoEpi":"MyoEpi"},"Hemolymph":{"1":"GH-like (c1)","6":"MH1-like (c6)","11":"cycling/PH-like (c11)","10":"c10 (unknown)","9":"c9 (unknown)"},"Hepatopancreas":{"6":"proliferative (c6)","2":"epithelial/immune (c2)","7":"R-like (c7)","15":"neuroendocrine (c15)","14":"c14 (unknown)"}}

def bezier(x0,y0,x1,y1):
 t=np.linspace(0,1,80); x=(1-t)**3*x0+3*(1-t)**2*t*(x0+0.22*(x1-x0))+3*(1-t)*t**2*(x1-0.22*(x1-x0))+t**3*x1
 y=(1-t)**3*y0+3*(1-t)**2*t*y0+3*(1-t)*t**2*y1+t**3*y1
 return x,y

def main():
 # only supported/trend links
 n_total=len(P)
 S=P.dropna(subset=["score","p_perm"]).copy(); n_valid=len(S)
 S=S[S.q_perm<0.10].copy(); n_supported=len(S)
 print(f"permutation rows: total={n_total}, valid={n_valid}, supported/trend={n_supported}")
 progs=[p for p in ["Progenitor","Macrophage-like","Hemocytes","ECM-producing","Capsular","Epithelial","SLC-rich","Lymphocyte precursor-like"] if p in set(S.program)]
 # programs without supported links omitted entirely
 left_y={p:0.78-0.20*i for i,p in enumerate(progs)}
 mid_y={"Hepatic sinus":0.78,"Hemolymph":0.50,"Hepatopancreas":0.22}
 # unique right labels from supported links
 right=[]; labels={}
 right_p={}
 for o in ORGS:
  sub=S[S.organ==o].copy(); sub["label"]=sub.best_group.astype(str).map(CELL[o]).fillna(sub.best_group.astype(str)+" (unknown)")
  labs=list(dict.fromkeys(sub["label"].tolist())); labels[o]=labs
  ys=np.linspace(mid_y[o]+0.07,mid_y[o]-0.07,len(labs)) if len(labs)>1 else [mid_y[o]]
  for l,y in zip(labs,ys):
   right.append((o,l,y))
   right_p[(o,l)]=float(sub[sub.label==l].p_perm.min())
 rpos={(o,l):y for o,l,y in right}
 fig=plt.figure(figsize=(5.7,4.0),facecolor="white"); ax=fig.add_axes([0.02,0.11,0.96,0.79]); ax.axis("off"); ax.set_xlim(0,1); ax.set_ylim(0,1)
 # nodes
 for p,y in left_y.items():
  ax.scatter([0.075],[y],s=24,c="#555555",zorder=8); ax.text(0.06,y,p,ha="right",va="center",fontsize=5.2)
 for o,y in mid_y.items():
  ax.scatter([0.30],[y],s=42,c=COL[o],zorder=8); ax.text(0.30,y+0.055,O2[ORGS.index(o)],ha="center",va="bottom",fontsize=5.0,color=COL[o],weight="bold")
 for o,l,y in right:
  ax.scatter([0.51],[y],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8)
  ax.text(0.525,y,"%s  p=%.3g"%(l,right_p[(o,l)]),ha="left",va="center",fontsize=5.0,color=COL[o])
 # draw each supported flow with identical style on both segments
 for p in progs:
  sub=S[S.program==p].sort_values("score",ascending=False)
  n=len(sub)
  for k,r in enumerate(sub.itertuples()):
   q=float(r.q_perm); lw=2.2+1.5*max(0,float(r.score))/4 if q<0.05 else 1.6
   ls="-" if q<0.05 else (0,(4,2)); alpha=0.9 if q<0.05 else 0.75
   # left segment
   x,y=bezier(0.075,left_y[p],0.30,mid_y[r.organ])
   ax.plot(x,y,lw=lw,ls=ls,color=COL[r.organ],alpha=alpha,zorder=2,solid_capstyle="round")
   # right segment to specific cell label
   lab=str(r.best_group); lab=CELL[r.organ].get(lab,lab+" (unknown)")
   yy=rpos[(r.organ,lab)]
   x,y=bezier(0.30,mid_y[r.organ],0.51,yy)
   ax.plot(x,y,lw=lw,ls=ls,color=COL[r.organ],alpha=alpha,zorder=2,solid_capstyle="round")

 ax.set_title("Statistically supported LyO programme sharing",fontsize=6.5,pad=4)
 fig.legend(handles=[Patch(color=COL[o],label=O2[i]) for i,o in enumerate(ORGS)],loc="lower center",bbox_to_anchor=(0.5,0.012),ncol=3,frameon=False,fontsize=5.0,handlelength=0.8,columnspacing=0.7,title="Organ",title_fontsize=5.0)
 fig.legend(handles=[Line2D([0],[0],color="#333333",lw=2.2,label="FDR < 0.05"),Line2D([0],[0],color="#333333",lw=1.5,ls=(0,(4,2)),label="trend 0.05–0.10")],loc="lower left",bbox_to_anchor=(0.03,0.012),ncol=1,frameon=False,fontsize=5.0,handlelength=1.5,labelspacing=0.25)
 fig.text(0.50,0.075,"Only supported/trend links shown; unsupported programmes omitted. Permutation p values are shown beside the right-hand cell types. Marker-level sharing, not migration.",ha="center",fontsize=5.0,color="#555555")
 p=os.path.join(PREV,"Fig6_alluvial_supported_only_rightp_close_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); print("saved",p)
 pqa=os.path.join(PREV,"Fig6_alluvial_supported_only_rightp_close_qa.pdf"); fig.savefig(pqa,bbox_inches="tight",pad_inches=0.04,facecolor="white"); print("saved",pqa)
 deliv=os.path.join(r"D:\桌面\single cell\文章配图","最终_2026-10-19","四器官模型_河流图"); os.makedirs(deliv,exist_ok=True)
 stem=os.path.join(deliv,"Fig6_alluvial_supported_only_rightp")
 fig.savefig(stem+".svg",bbox_inches="tight",pad_inches=0.04,facecolor="white")
 fig.savefig(stem+".tiff",dpi=1200,bbox_inches="tight",pad_inches=0.04,facecolor="white",pil_kwargs={"compression":"tiff_lzw"})
 with Image.open(stem+".tiff") as im:
  if "A" in im.getbands():
   rgba=im.convert("RGBA"); bg=Image.new("RGB",rgba.size,(255,255,255)); bg.paste(rgba,mask=rgba.getchannel("A")); tmp=stem+".tmp.tiff"; bg.save(tmp,format="TIFF",dpi=(1200,1200),compression="tiff_lzw"); os.replace(tmp,stem+".tiff")
 print("saved",stem+".svg"); print("saved",stem+".tiff")
 plt.close(fig)

if __name__=="__main__": main()
