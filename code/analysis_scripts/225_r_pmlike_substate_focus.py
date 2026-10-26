# -*- coding: utf-8 -*-
"""225_r_pmlike_substate_focus.py"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys, numpy as np, pandas as pd
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
def main():
 d=pd.read_csv(os.path.join(OUT,'Fig5_dynamic_subpop_subcluster_composition.csv'))
 d['key']=d.population.astype(str)+'_s'+d.subcluster.astype(str)
 keys=['R_s0','R_s2','R_s3','PM-like_s0','PM-like_s2']
 d=d[d.key.isin(keys)].copy(); d['order']=d.key.map({k:i for i,k in enumerate(keys)}); d=d.sort_values('order')
 fig=plt.figure(figsize=(4.9,3.6),facecolor='white'); ax=fig.add_axes([0.30,0.23,0.58,0.66])
 y=np.arange(len(d))[::-1]
 for yy,r in zip(y,d.itertuples()):
  ax.plot([r.pct_control,r.pct_hemospasia],[yy,yy],color='#BDBDBD',lw=1.0,zorder=1)
  ax.scatter([r.pct_control],[yy],s=20,c='#4C78A8',edgecolors='white',linewidths=.35,zorder=3)
  ax.scatter([r.pct_hemospasia],[yy],s=20,c='#D6604D',edgecolors='white',linewidths=.35,zorder=3)
  ax.text(103,yy,f"{r.fold_hemospasia_vs_control:.2f}×",ha='left',va='center',fontsize=5.1,color='#B2182B' if r.fold_hemospasia_vs_control>1 else '#2166AC')
 ax.set_yticks(y); ax.set_yticklabels([f"{r.key} ({r.n_control}/{r.n_hemospasia})" for r in d.itertuples()],fontsize=5.2)
 ax.set_xlim(0,116); ax.set_ylim(-.55,len(d)-.45); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel('Percentage within parent state (%)',fontsize=5.5)
 ax.tick_params(axis='x',labelsize=5.0,length=1.4,width=.4); ax.tick_params(axis='y',length=0,pad=2)
 for sp in ('top','right'): ax.spines[sp].set_visible(False)
 ax.grid(axis='x',color='#EEEEEE',lw=.35,zorder=0)
 fig.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#4C78A8',markersize=4,label='Control'),Line2D([0],[0],marker='o',color='none',markerfacecolor='#D6604D',markersize=4,label='Hemospasia')],loc='lower center',bbox_to_anchor=(.5,.085),ncol=2,frameon=False,fontsize=5.0,handlelength=.7,columnspacing=1.0)
 ax.set_title('R and PM-like substate composition',fontsize=6.8,pad=5)
 fig.text(.5,.018,'Fold change at right. PM-like total cells: 583 → 127; only substates with >=20 cells and >=5 per group are shown.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_R_PMlike_substate_focus_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_R_PMlike_substate_focus.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_R_PMlike_substate_focus.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
