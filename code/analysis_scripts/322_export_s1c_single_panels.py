# -*- coding: utf-8 -*-
"""Export S1-C as a single UCell panel and retain mapping/overlap as separate panels."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle, Patch
from PIL import Image,ImageDraw,ImageFont

ROOT=Path(r'D:\桌面\single cell'); OUT=ROOT/'analysis_output'; FINAL=OUT/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_final'; D=FINAL/'S1C_individual'; D.mkdir(parents=True,exist_ok=True)
assert fm.findfont('Times New Roman',fallback_to_default=False)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'font.size':8,'axes.titlesize':9,'axes.labelsize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':0.65,'savefig.facecolor':'white','svg.fonttype':'none','pdf.fonttype':42})
def save(fig,stem):
 fig.savefig(stem.with_suffix('.png'),dpi=300,bbox_inches='tight'); fig.savefig(stem.with_suffix('.svg'),format='svg',bbox_inches='tight'); fig.savefig(stem.with_suffix('.pdf'),format='pdf',bbox_inches='tight'); fig.savefig(stem.with_suffix('.tiff'),format='tiff',dpi=1200,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); plt.close(fig); print(stem)

# S1-C canonical: LyO progenitor UCell lollipop.
u=pd.read_csv(OUT/'figures_pub'/'Supplementary_niche'/'Supp14_LyO_progenitor_Hs_UCell.csv')
order=['PH','GH','TS','MH2','F','R','B','M','MyoEpi','PM-like']; u=u.set_index(u.iloc[:,0]).reindex(order)
mean=u['mean'].astype(float).values; frac=u['frac_above_global'].astype(float).values
fig,ax=plt.subplots(figsize=(5.2,3.8)); y=np.arange(len(order)); colors=['#B2182B' if x=='PH' else '#8A8A8A' for x in order]
ax.hlines(y,0.46,mean,color='#D5D5D5',lw=1.1); ax.scatter(mean,y,s=22+frac*70,c=colors,zorder=3,edgecolors='#333',linewidths=.35)
for i,(m,f) in enumerate(zip(mean,frac)): ax.text(m+.003,i,f'{m:.3f} / {100*f:.0f}%',va='center',fontsize=6.5,color='#B2182B' if order[i]=='PH' else '#444')
ax.set_yticks(y); ax.set_yticklabels(order); ax.invert_yaxis(); ax.set_xlim(.46,.55); ax.set_xlabel('Mean LyO progenitor UCell score'); ax.set_title('LyO progenitor program in Hs',pad=5)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,D/'S1_C_LyO_progenitor_UCell')

# Extra single panel: LyO cluster-to-Hs marker map.
heat=pd.read_csv(OUT/'lyo_cross_organ_mapping_scores.csv'); heat=heat[heat.organ=='Hepatic sinus']; best=pd.read_csv(OUT/'Fig_lyo_cluster_best_matches.csv'); nm=dict(zip(best.lyo_cluster,best.lyo_name)); cell_order=['PH','GH','TS','MH2','F','R','B','M','MyoEpi','PM-like']; mat=heat.pivot(index='lyo_cluster',columns='target_group',values='z_score').reindex(columns=cell_order).sort_index(); labs=[f'{int(c)} {nm.get(int(c),"")}' for c in mat.index]
fig,ax=plt.subplots(figsize=(5.7,5.6)); vmax=np.nanmax(np.abs(mat.values)); im=ax.imshow(mat.values,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax); ax.set_xticks(range(len(cell_order))); ax.set_xticklabels(cell_order,rotation=55,ha='right'); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=6.5); ax.set_xlabel('Hs cell state'); ax.set_ylabel('Lymphoid-organ cluster'); ax.set_title('LyO cluster-to-Hs marker mapping',pad=5)
if 13 in mat.index:
 r=list(mat.index).index(13); ax.add_patch(Rectangle((-.5,r-.5),len(cell_order),1,fill=False,ec='#111',lw=1))
for sp in ax.spines.values(): sp.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03); cb.set_label('Marker-overlap z',fontsize=7); cb.ax.tick_params(labelsize=6)
fig.tight_layout(); save(fig,D/'S1_C_extra_LyO_cluster_mapping')

# Extra single panel: marker-level overlap.
d=pd.read_csv(OUT/'figures_pub'/'Supplementary_niche'/'Supp14_LyO_progenitor_overlap.csv'); d['label']=d['label'].replace({'PH':'PH','Niche':'Niche'}); oo=['PH','Niche','TS','MH2','MyoEpi','M']; d=d.set_index('label').reindex(oo).reset_index(); x=range(len(d))
fig,ax=plt.subplots(figsize=(4.6,4.0)); ax.scatter(x,d['fold'],s=[75 if v=='PH' else 24 for v in d.label],c=['#B2182B' if v=='PH' else '#8A8A8A' for v in d.label],zorder=3,edgecolors='#333',linewidths=.4)
for i,r in d.iterrows():
 if r['fold']>0: ax.annotate(f"{int(r['overlap'])} genes; {r['fold']:.1f}×\nP={r['p_hypergeom']:.1e}",(i,r['fold']),xytext=(8,4),textcoords='offset points',fontsize=6.5,color='#B2182B',ha='left',va='bottom')
 else: ax.text(i,1.5,'0 genes\nP=1',ha='center',va='bottom',fontsize=6.2,color='#666')
ax.set_xticks(list(x)); ax.set_xticklabels(d.label,rotation=40,ha='right'); ax.set_ylim(0,80); ax.set_xlim(-.6,5.6); ax.set_ylabel('Fold enrichment'); ax.set_title('LyO progenitor top50 overlap',pad=5)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,D/'S1_C_extra_marker_overlap')

# Recompose S1 with A, B and canonical UCell C.
def load(p):
 im=Image.open(p).convert('RGBA'); bg=Image.new('RGBA',im.size,'white'); bg.alpha_composite(im); return bg.convert('RGB')
def contain(im,w,h): im=im.copy(); im.thumbnail((w,h),Image.Resampling.LANCZOS); return im
A=load(FINAL/'S1_A_PH_substate_markers.png'); B=load(FINAL/'S1_B_PH_identity_vs_cellcycle.png'); C=load(D/'S1_C_LyO_progenitor_UCell.png'); W=2400; margin=40; gap=30; header=100; top_h=1000; colw=(W-2*margin-gap)//2; A=contain(A,colw,top_h); B=contain(B,colw,top_h); C=contain(C,1800,1200); y1=margin+header; y2=y1+top_h+gap; H=y2+C.height+margin
S=Image.new('RGB',(W,H),'white'); dr=ImageDraw.Draw(S); dr.text((W//2,35),'Supplementary Figure S1 | PH substate identity evidence',font=ImageFont.truetype(r'C:\Windows\Fonts\timesbd.ttf',38),fill='#0f172a',anchor='mm'); S.paste(A,(margin+(colw-A.width)//2,y1)); S.paste(B,(margin+colw+gap+(colw-B.width)//2,y1)); S.paste(C,((W-C.width)//2,y2))
stem=FINAL/'Supplementary_Figure_S1'; S.save(stem.with_suffix('.png'),dpi=(300,300)); S.save(stem.with_suffix('.tiff'),dpi=(1200,1200),compression='tiff_lzw'); S.save(stem.with_suffix('.pdf'),'PDF',resolution=300,dpi=(300,300)); print(stem)
