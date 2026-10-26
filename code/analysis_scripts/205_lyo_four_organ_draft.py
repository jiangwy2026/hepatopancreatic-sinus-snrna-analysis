# -*- coding: utf-8 -*-
"""205_lyo_four_organ_draft.py
Draft second figure: LyO clusters/markers + cross-organ programmes + model.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
LYO=r"D:\桌面\single cell\DeGene.list.xls"; SC=os.path.join(OUT,"lyo_cross_organ_mapping_scores.csv")
NAME={0:"ECM-producing",1:"Macrophage-like 1",2:"Macrophage-like 2",3:"Unknown",4:"Lymphocyte precursor-like",5:"Capsular",6:"Unknown",7:"Hemocytes",8:"SLC-rich",9:"Unknown",10:"Unknown",11:"Epithelial",12:"Unknown",13:"Progenitor",14:"Unknown",15:"Unknown"}
PROGRAMS=[13,2,7,0,5,11,8,4]
PNAME={13:"Progenitor",2:"Macrophage-like",7:"Hemocytes",0:"ECM-producing",5:"Capsular",11:"Epithelial",8:"SLC-rich",4:"Lymphocyte precursor-like"}
TARGETS=[("Hepatic sinus",["PH","MH2","GH","R","PM-like","MyoEpi"]),("Hemolymph",["11","6","1","10","9"]),("Hepatopancreas",["6","2","7","14","15"])]

def main():
 d=pd.read_excel(LYO,sheet_name="Top50.DeGene.list"); d["Target_Cluster"]=d.Target_Cluster.astype(int)
 top={}
 for c,g in d.groupby("Target_Cluster"):
  g=g.sort_values("Log2FC",ascending=False)
  top[c]=g.head(5).copy()
 # selected marker table
 rows=[]
 for c in range(16):
  for r in top[c].itertuples(): rows.append({"cluster":c,"name":NAME[c],"gene":r._3,"log2FC":float(r.Log2FC)})
 pd.DataFrame(rows).to_csv(os.path.join(OUT,"Fig6_LyO_top_markers.csv"),index=False,encoding="utf-8-sig")
 # cross-organ matrix
 S=pd.read_csv(SC)
 cols=[]; labels=[]; groups=[]
 for org,gs in TARGETS:
  for g in gs:
   x=S[(S.organ==org)&(S.target_group==g)]
   if len(x): cols.append(x.iloc[0]["z_score"] if len(x)==1 else np.nan)
   labels.append(g); groups.append(org)
 # build exact matrix from program rows
 M=[]
 for c in PROGRAMS:
  row=[]
  for org,g in zip(groups,labels):
   x=S[(S.organ==org)&(S.target_group==g)&(S.lyo_cluster==c)]
   row.append(float(x.iloc[0].z_score) if len(x) else np.nan)
  M.append(row)
 M=np.array(M); pd.DataFrame(M,index=[PNAME[c] for c in PROGRAMS],columns=["%s:%s"%(o,g) for o,g in zip(groups,labels)]).to_csv(os.path.join(OUT,"Fig6_cross_organ_program_matrix.csv"),encoding="utf-8-sig")
 fig=plt.figure(figsize=(7.2,6.4),facecolor="white")
 # panel a table
 axa=fig.add_axes([0.02,0.39,0.55,0.55]); axa.axis("off")
 axa.set_title("Lymphoid organ clusters and top markers",fontsize=6.7,pad=4)
 axa.text(0.0,0.98,"Cluster / annotation",fontsize=5.0,weight="bold",va="top")
 axa.text(0.39,0.98,"Top 5 markers by log2FC",fontsize=5.0,weight="bold",va="top")
 for i,c in enumerate(range(16)):
  y=0.91-i*0.057
  axa.text(0.0,y,"%d  %s"%(c,NAME[c]),fontsize=4.8,va="center",color="#B2182B" if c in (0,1,2,4,5,7,8,11,13) else "#666666")
  genes=" · ".join(top[c]["Gene Name"].astype(str).tolist())
  axa.text(0.39,y,genes,fontsize=4.6,va="center",wrap=True)
 axa.set_xlim(0,1); axa.set_ylim(0,1)
 # panel b heatmap
 axb=fig.add_axes([0.63,0.44,0.34,0.49])
 im=axb.imshow(M,aspect="auto",cmap="RdBu_r",vmin=-3,vmax=3)
 axb.set_yticks(range(len(PROGRAMS))); axb.set_yticklabels([PNAME[c] for c in PROGRAMS],fontsize=4.8)
 axb.set_xticks(range(len(labels))); axb.set_xticklabels(labels,fontsize=4.3,rotation=45,ha="right")
 axb.tick_params(length=0)
 for i in range(len(PROGRAMS)):
  for j in range(len(labels)):
   v=M[i,j]
   if np.isfinite(v): axb.text(j,i,"%.1f"%v,ha="center",va="center",fontsize=3.9,color="white" if abs(v)>1.5 else "#222222")
 # organ labels
 for j0,j1,org in [(0,5,"Hepatic sinus"),(6,10,"Hemolymph"),(11,15,"Hepatopancreas")]:
  axb.text((j0+j1)/2,-1.7,org,ha="center",va="bottom",fontsize=4.7,weight="bold")
 axb.set_title("Cross-organ programme score",fontsize=6.5,pad=20)
 for sp in axb.spines.values(): sp.set_visible(False)
 cb=fig.colorbar(im,ax=axb,fraction=0.030,pad=0.02); cb.set_label("signature z-score",fontsize=4.3); cb.ax.tick_params(labelsize=4.0,length=0.6); cb.outline.set_visible(False)
 # panel c model
 axc=fig.add_axes([0.04,0.02,0.92,0.28]); axc.set_xlim(0,1); axc.set_ylim(0,1); axc.axis("off")
 def box(x,y,w,h,text,fc,ec,fc_text="white"):
  axc.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012,rounding_size=0.02",lw=0.7,ec=ec,fc=fc,zorder=2))
  axc.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=4.9,color=fc_text,zorder=3)
 def arr(p1,p2,style="-|>",ls="-",col="#333333"):
  axc.add_patch(FancyArrowPatch(p1,p2,arrowstyle=style,mutation_scale=8,lw=1.0,color=col,ls=ls,zorder=4,shrinkA=2,shrinkB=2))
 box(0.02,0.34,0.17,0.30,"Hepatopancreas","#F4E1C1","#8C6D3A","#4A3418")
 box(0.26,0.34,0.17,0.30,"Hepatic sinus","#F2C6C2","#B2182B","#5A1111")
 box(0.55,0.60,0.17,0.26,"Lymphoid organ","#D9E8F5","#2166AC","#123A5C")
 box(0.78,0.34,0.17,0.30,"Hemolymph","#D9E8F5","#2166AC","#123A5C")
 arr((0.19,0.49),(0.26,0.49),ls=(0,(3,2)),col="#8C6D3A")
 axc.text(0.225,0.55,"ECM / metabolic\nsignals",ha="center",fontsize=4.2,color="#8C6D3A")
 arr((0.43,0.54),(0.55,0.66),ls=(0,(3,2)),col="#2166AC")
 arr((0.55,0.64),(0.43,0.52),ls=(0,(3,2)),col="#2166AC")
 axc.text(0.49,0.69,"progenitor / hemocyte\nprogrammes",ha="center",fontsize=4.1,color="#2166AC")
 arr((0.43,0.42),(0.78,0.42),col="#B2182B")
 axc.text(0.61,0.35,"PH / GH / MH2 output",ha="center",fontsize=4.2,color="#B2182B")
 arr((0.72,0.60),(0.78,0.50),ls=(0,(3,2)),col="#2166AC")
 axc.text(0.81,0.61,"immune\nmaturation",ha="center",fontsize=4.1,color="#2166AC")
 axc.text(0.02,0.08,"Solid: supported by known anatomy/function. Dashed: marker-level candidate relationship. HPT not included (insufficient direct evidence in L. vannamei).",fontsize=4.4,color="#555555")
 fig.text(0.5,0.985,"Four-organ relationship and candidate hematopoietic model",ha="center",va="top",fontsize=7.0)
 p=os.path.join(PREV,"Fig6_four_organ_draft_preview.png"); fig.savefig(p,dpi=350,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)
if __name__=="__main__": main()
