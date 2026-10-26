# -*- coding: utf-8 -*-
"""231_group_specific_focused_lr.py"""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
sc.settings.verbosity=0
G={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
SENDERS=['F','R','B','M','MyoEpi','PM-like']; RECEIVERS=['PH','GH','TS','MH2']
PAIRS=[('Wnt5b','ncbi_113806005','Fzd4','ncbi_113813183','Wnt'),('Wnt6','ncbi_113812192','Fzd1','ncbi_113813182','Wnt'),('PDGF/VEGF-like','ncbi_113804544','Flt1','ncbi_113803634','VEGF/PDGF'),('FGF16','ncbi_113829301','FGFR1','ncbi_113803204','FGF'),('THBS4','ncbi_113824160','ITGA5','ncbi_113807444','ECM-integrin'),('THBS4','ncbi_113824160','ITGB3','ncbi_113805718','ECM-integrin'),('THBS4','ncbi_113824160','SDC2','ncbi_113807065','ECM-syndecan'),('Igfbp7','ncbi_113811702','Igf1r','ncbi_113817108','IGF')]
def main():
 A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad')); A.X=A.layers['counts'].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 A.obs['ct']=A.obs.vendor_cluster.astype(str).map(G); A.obs['grp']=A.obs['sample'].astype(str).map({'blank':'Control','treated':'Hemospasia'}); gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
 rows=[]
 for lname,lgid,rname,rgid,path in PAIRS:
  if lgid not in gid2i or rgid not in gid2i: continue
  lv=A[:,gid2i[lgid]].X.toarray().ravel(); rv=A[:,gid2i[rgid]].X.toarray().ravel()
  for grp in ['Control','Hemospasia']:
   gm=A.obs.grp.astype(str).values==grp
   for s in SENDERS:
    ms=gm&(A.obs.ct.astype(str).values==s); lm=float(lv[ms].mean()) if ms.any() else 0; lp=100*float((lv[ms]>0).mean()) if ms.any() else 0
    for r in RECEIVERS:
     mr=gm&(A.obs.ct.astype(str).values==r); rm=float(rv[mr].mean()) if mr.any() else 0; rp=100*float((rv[mr]>0).mean()) if mr.any() else 0
     score=np.sqrt(lm*rm) if lp>=5 and rp>=5 else np.nan
     rows.append({'pathway':path,'ligand':lname,'receptor':rname,'sender':s,'receiver':r,'group':grp,'ligand_pct':lp,'receptor_pct':rp,'score':score})
 R=pd.DataFrame(rows); R.to_csv(os.path.join(OUT,'Fig5_focused_LR_by_group.csv'),index=False,encoding='utf-8-sig')
 C=R.pivot_table(index=['pathway','ligand','receptor','sender','receiver'],columns='group',values='score',aggfunc='first').reset_index()
 C['delta']=C['Hemospasia']-C['Control']; C=C.dropna(subset=['Control','Hemospasia']); C['absdelta']=C.delta.abs(); C=C.sort_values('absdelta',ascending=False).head(20).sort_values('delta')
 C.to_csv(os.path.join(OUT,'Fig5_focused_LR_delta_top.csv'),index=False,encoding='utf-8-sig')
 fig=plt.figure(figsize=(5.6,4.5),facecolor='white'); ax=fig.add_axes([0.47,0.20,0.49,0.70]); y=np.arange(len(C))
 ax.hlines(y,C.Control,C.Hemospasia,color='#BBBBBB',lw=1.0,zorder=1)
 ax.scatter(C.Control,y,s=16,c='#4C78A8',edgecolors='white',linewidths=.35,zorder=3,label='Control')
 ax.scatter(C.Hemospasia,y,s=16,c='#D6604D',edgecolors='white',linewidths=.35,zorder=3,label='Hemospasia')
 labels=[f"{r.ligand} {r.sender} → {r.receptor} {r.receiver}" for r in C.itertuples()]
 ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=4.7); ax.set_xlabel('Candidate interaction score',fontsize=5.3); ax.set_xlim(left=0)
 ax.tick_params(axis='x',labelsize=4.9,length=1.2); ax.tick_params(axis='y',length=0,pad=2)
 for sp in ('top','right'): ax.spines[sp].set_visible(False)
 ax.grid(axis='x',color='#EEEEEE',lw=.35); ax.set_axisbelow(True)
 ax.set_title('Group-specific focused ligand–receptor screen',fontsize=6.8,pad=5)
 fig.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#4C78A8',markersize=4,label='Control'),Line2D([0],[0],marker='o',color='none',markerfacecolor='#D6604D',markersize=4,label='Hemospasia')],loc='lower center',bbox_to_anchor=(.5,.075),ncol=2,frameon=False,fontsize=5.0)
 fig.text(.5,.015,'Orthology-based candidate screen; not a validated shrimp ligand–receptor database. Top 20 interactions by absolute change.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_focused_LR_group_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_focused_LR_group.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_focused_LR_group.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
