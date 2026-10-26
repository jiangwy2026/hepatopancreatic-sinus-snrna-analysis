# -*- coding: utf-8 -*-
"""228_deg_sharing_upset.py
UpSet-style summary of shared and unique DEGs across dynamic populations."""
import os, sys, itertools
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化'); DE=os.path.join(OUT,'dynamic_subpop_de')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
SETS=['PH','MH2','R','PM-like']; COL={'PH':'#B2182B','MH2':'#8C510A','R':'#41AB5D','PM-like':'#7B3294'}
def main():
 sig={}
 for s in SETS:
  d=pd.read_csv(os.path.join(DE,f'DE_{s}_treated_vs_control.csv'))
  d=d[(d.pvals_adj<0.05)&(d.logfoldchanges.abs()>0.5)]
  sig[s]=set(d.sym.dropna().astype(str))
 combos=[]; rows=[]
 for k in range(1,5):
  for comb in itertools.combinations(SETS,k):
   exact=set.intersection(*[sig[s] for s in comb])
   for s in SETS:
    if s not in comb: exact-=sig[s]
   combos.append((comb,len(exact)))
   for g in sorted(exact): rows.append({'combination':'+'.join(comb),'gene':g,'n_sets':k})
 combos=sorted(combos,key=lambda x:x[1],reverse=True)[:18]
 pd.DataFrame(rows).to_csv(os.path.join(OUT,'Fig5_DEG_shared_unique.csv'),index=False,encoding='utf-8-sig')
 fig=plt.figure(figsize=(7.2,4.2),facecolor='white')
 axb=fig.add_axes([0.18,0.62,0.76,0.25]); axm=fig.add_axes([0.18,0.18,0.76,0.38]); axs=fig.add_axes([0.04,0.18,0.05,0.38])
 n=len(combos); xs=np.arange(n); sizes=[c[1] for c in combos]
 axb.bar(xs,sizes,color=['#444444' if c[1]>0 else '#CCCCCC' for c in combos],width=.72)
 for x,v in zip(xs,sizes): axb.text(x,v+max(sizes)*.02,str(v),ha='center',va='bottom',fontsize=4.6)
 axb.set_xlim(-.6,n-.4); axb.set_ylim(0,max(sizes)*1.18); axb.set_xticks([]); axb.set_ylabel('Intersection size',fontsize=5.2); axb.tick_params(axis='y',labelsize=5.0,length=1.2)
 for sp in ('top','right'): axb.spines[sp].set_visible(False)
 for j,(comb,size) in enumerate(combos):
  for i,s in enumerate(SETS):
   on=s in comb
   axm.scatter(j,len(SETS)-1-i,s=18 if on else 7,c=COL[s] if on else '#DDDDDD',edgecolors='white',linewidths=.3)
  if len(comb)>1:
   ys=[len(SETS)-1-SETS.index(s) for s in comb]; axm.plot([j,j],[min(ys),max(ys)],color='#777777',lw=.8,zorder=0)
 axm.set_xlim(-.6,n-.4); axm.set_ylim(-.5,len(SETS)-.5); axm.set_xticks(xs); axm.set_xticklabels([]); axm.set_yticks(range(len(SETS))); axm.set_yticklabels(['PH','MH2','R','PM'][::-1],fontsize=5.2)
 axm.tick_params(axis='y',length=0,pad=2); axm.tick_params(axis='x',length=0)
 for sp in axm.spines.values(): sp.set_visible(False)
 axs.barh(range(len(SETS)),[len(sig[s]) for s in SETS[::-1]],color=[COL[s] for s in SETS[::-1]],height=.55)
 axs.set_yticks([]); axs.set_xlabel('Set size',fontsize=5.0); axs.tick_params(axis='x',labelsize=4.7,length=1.0); axs.set_ylim(-.5,len(SETS)-.5)
 for sp in ('top','right'): axs.spines[sp].set_visible(False)
 fig.text(.5,.96,'Shared and unique DEGs across dynamic populations',ha='center',va='top',fontsize=7.2)
 fig.text(.5,.03,'DEGs: adjusted P < 0.05 and |log2FC| > 0.5; exact intersections shown. MyoEpi excluded by cell-number threshold.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_DEG_sharing_upset_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_DEG_sharing_upset.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_DEG_sharing_upset.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); print(pd.DataFrame([{'combination':'+'.join(c),'n':nv} for c,nv in combos]).to_string(index=False)); plt.close(fig)
if __name__=='__main__': main()
