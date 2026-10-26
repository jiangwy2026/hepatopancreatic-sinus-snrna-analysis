# -*- coding: utf-8 -*-
"""232_tf_expression_shift.py
TF expression shifts across dynamic populations (proxy; not SCENIC activity)."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化'); DE=os.path.join(OUT,'dynamic_subpop_de')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
PARENTS=['PH','MH2','R','PM-like']
def main():
 cand=pd.read_csv(os.path.join(OUT,'PH_TF_candidates.csv')).sort_values('evidence_score',ascending=False)
 tfs=cand.TF.astype(str).head(14).tolist()
 for x in ['Stat5b','pnt','foxo','br','Trl']:
  if x not in tfs: tfs.append(x)
 mats={}; qs={}
 for par in PARENTS:
  d=pd.read_csv(os.path.join(DE,f'DE_{par}_treated_vs_control.csv'))
  dm=d.drop_duplicates('sym').set_index('sym')
  mats[par]=[float(dm.loc[t,'logfoldchanges']) if t in dm.index else np.nan for t in tfs]
  qs[par]=[float(dm.loc[t,'pvals_adj']) if t in dm.index else np.nan for t in tfs]
 M=pd.DataFrame(mats,index=tfs); Q=pd.DataFrame(qs,index=tfs); M.to_csv(os.path.join(OUT,'Fig5_TF_expression_log2fc.csv'),encoding='utf-8-sig'); Q.to_csv(os.path.join(OUT,'Fig5_TF_expression_fdr.csv'),encoding='utf-8-sig')
 fig=plt.figure(figsize=(5.2,4.0),facecolor='white'); ax=fig.add_axes([0.26,0.18,0.62,0.70])
 im=ax.imshow(np.clip(M.values,-3,3),aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3)
 ax.set_xticks(range(len(PARENTS))); ax.set_xticklabels(PARENTS,rotation=38,ha='right',rotation_mode='anchor',fontsize=5.2)
 ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs,fontsize=5.1); ax.tick_params(length=0,pad=2)
 for i in range(len(tfs)):
  for j,par in enumerate(PARENTS):
   if np.isfinite(Q.loc[tfs[i],par]) and Q.loc[tfs[i],par]<0.05: ax.text(j,i,'*',ha='center',va='center',fontsize=6,color='#222222')
 for sp in ax.spines.values(): sp.set_visible(False)
 cax=fig.add_axes([0.90,0.28,0.018,0.42]); cb=fig.colorbar(im,cax=cax); cb.set_label('log2FC (Hemospasia/Control)',fontsize=5.0,labelpad=4); cb.ax.tick_params(labelsize=4.9,length=1.2,width=.4); cb.outline.set_linewidth(.4)
 fig.text(.5,.96,'TF expression shifts across dynamic populations',ha='center',va='top',fontsize=6.8)
 fig.text(.5,.03,'TF expression proxy; SCENIC regulon activity not available due incomplete motif database. * FDR < 0.05.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_TF_expression_shift_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_TF_expression_shift.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_TF_expression_shift.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
