import os, sys, numpy as np, pandas as pd
os.environ.setdefault('MPLCONFIGDIR', r'D:\桌面\single cell\analysis_output\.mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
R=r'D:\桌面\single cell'; O=os.path.join(R,'analysis_output'); P=os.path.join(O,'figures_pub','fig5_parts'); D=os.path.join(R,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化'); os.makedirs(P,exist_ok=True); os.makedirs(D,exist_ok=True)
sys.path.insert(0,O)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Nimbus Roman','Liberation Serif','DejaVu Serif'],'svg.fonttype':'none','pdf.fonttype':42})
C={'G1':'#BFBFBF','S':'#4393C3','G2M':'#B2182B'}
d=pd.read_csv(os.path.join(O,'Fig5_dynamic_cellcycle_by_group.csv')); keep=['PH','MH2','R','PM-like']; d=d[d.population.isin(keep)].copy()
order=[]
for g in keep: order += [(g,'Control'),(g,'Hemospasia')]
y=[]; cur=0
for g in keep: y += [cur+.18,cur-.18]; cur-=1
fig=plt.figure(figsize=(4.7,3.2),facecolor='white'); ax=fig.add_axes([.34,.20,.58,.68])
for (g,grp),yy in zip(order,y):
 r=d[(d.population==g)&(d.group==grp)].iloc[0]; left=0
 for ph in ['G1','S','G2M']:
  v=float(r[f'pct_{ph}']); ax.barh(yy,v,left=left,color=C[ph],height=.34,edgecolor='white',linewidth=.35)
  if v>=10: ax.text(left+v/2,yy,f'{v:.0f}',ha='center',va='center',fontsize=4.6,color='#333' if ph=='G1' else 'white')
  left+=v
ax.set_yticks(y); ax.set_yticklabels([f'{g} {grp}' for g,grp in order],fontsize=5.0); ax.set_xlim(0,100); ax.set_xlabel('Cells per phase (%)',fontsize=5.3); ax.tick_params(axis='x',labelsize=4.9,length=1.2); ax.tick_params(axis='y',length=0,pad=2)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.grid(axis='x',color='#EEEEEE',lw=.35); ax.set_axisbelow(True); ax.legend(handles=[Patch(color=C[p],label=p) for p in ['G1','S','G2M']],loc='upper center',bbox_to_anchor=(.5,1.17),ncol=3,frameon=False,fontsize=5.0)
ax.set_title('Cell-cycle phase by group',fontsize=6.6,pad=5)
fig.text(.5,.025,'MyoEpi excluded by cell-number threshold.',ha='center',va='bottom',fontsize=5.0,color='#555')
stem=os.path.join(P,'Fig5_cellcycle_dynamic_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
fig.savefig(os.path.join(D,'Fig5_cellcycle_dynamic.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(D,'Fig5_cellcycle_dynamic.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
print('saved',stem); plt.close(fig)
