# -*- coding: utf-8 -*-
"""Compact preview for LyO progenitor marker overlap."""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
ROOT=Path(r'D:\桌面\single cell'); OUT=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_previews'
OUT.mkdir(parents=True,exist_ok=True)
assert fm.findfont('Times New Roman',fallback_to_default=False)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'font.size':8,'axes.titlesize':9,'axes.labelsize':9.5,'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':0.7,'svg.fonttype':'none','savefig.facecolor':'white'})
d=pd.read_csv(ROOT/'analysis_output'/'figures_pub'/'Supplementary_niche'/'Supp14_LyO_progenitor_overlap.csv')
d['label']=d['label'].replace({'PH':'PH','Niche':'Niche'})
order=['PH','Niche','TS','MH2','MyoEpi','M']; d=d.set_index('label').reindex(order).reset_index(); x=range(len(d))
fig,ax=plt.subplots(figsize=(4.5,3.9)); cols=['#B2182B' if v=='PH' else '#8A8A8A' for v in d.label]
ax.scatter(x,d['fold'],s=[75 if v=='PH' else 24 for v in d.label],c=cols,zorder=3,edgecolors='#333',linewidths=0.4)
for i,r in d.iterrows():
    if r['fold']>0:
        ax.annotate(f"{int(r['overlap'])} genes; {r['fold']:.1f}×\nP={r['p_hypergeom']:.1e}",(i,r['fold']),xytext=(8,4),textcoords='offset points',fontsize=6.5,color='#B2182B',ha='left',va='bottom')
    else:
        ax.text(i,1.5,'0 genes\nP=1',ha='center',va='bottom',fontsize=6.2,color='#666666')
ax.set_xticks(list(x)); ax.set_xticklabels(d.label,rotation=40,ha='right'); ax.set_ylim(0,80); ax.set_xlim(-0.6,5.6)
ax.set_ylabel('Fold enrichment'); ax.set_title('LyO progenitor top50 overlap',pad=5)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); p=OUT/'Supp_LyO_overlap_compact_preview.png'; fig.savefig(p,dpi=300,bbox_inches='tight'); plt.close(fig); print(p)
