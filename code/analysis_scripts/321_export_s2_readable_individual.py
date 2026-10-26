# -*- coding: utf-8 -*-
"""Export individual S2 panels using the readability-optimized style."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

ROOT=Path(r'D:\桌面\single cell'); OUT=ROOT/'analysis_output'; STEP1=OUT/'step1'; STEP3=OUT/'step3'
DEST=OUT/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_final'/'readable_single_panels'; DEST.mkdir(parents=True,exist_ok=True)
assert fm.findfont('Times New Roman',fallback_to_default=False)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'font.size':8,'axes.titlesize':9,'axes.labelsize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':0.65,'savefig.facecolor':'white','svg.fonttype':'none','pdf.fonttype':42})
lr=pd.read_csv(STEP1/'step1_lr_communication_sonicparanoid.csv'); gdf=pd.read_csv(OUT/'figures_pub'/'Supplementary_niche'/'Supp7_group_LR.csv'); a=ad.read_h5ad(STEP3/'step3_lyo_scores.h5ad'); T=a.obsm['X_tsne']
g2v={str(g).replace('ncbi_',''):g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers['lognorm'].todense()).mean(1)
itgb1=expr(['113830322','113805715'])
def save(fig,stem):
 fig.savefig(stem.with_suffix('.png'),dpi=300,bbox_inches='tight'); fig.savefig(stem.with_suffix('.svg'),format='svg',bbox_inches='tight'); fig.savefig(stem.with_suffix('.pdf'),format='pdf',bbox_inches='tight'); fig.savefig(stem.with_suffix('.tiff'),format='tiff',dpi=1200,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); plt.close(fig); print(stem)

# A: faceted ECM-integrin / Notch heatmaps
src=['TS','MH2','MyoEpi','M']; tgt=['PH','GH']; fig=plt.figure(figsize=(5.5,3.8)); gs=fig.add_gridspec(1,2,wspace=0.5,left=0.12,right=0.95,bottom=0.20,top=0.84)
for k,axis in enumerate(['ECM-integrin','Notch']):
 ax=fig.add_subplot(gs[0,k]); piv=lr[lr.axis==axis].groupby(['source','target'])['prob'].sum().unstack(fill_value=0).reindex(index=src,columns=tgt,fill_value=0); im=ax.imshow(piv.values,cmap='Reds' if k==0 else 'Blues',aspect='auto'); ax.set_xticks([0,1]); ax.set_xticklabels(tgt); ax.set_yticks(range(4)); ax.set_yticklabels(src); ax.tick_params(length=0); ax.set_title(axis)
 for i in range(4):
  for j in range(2): ax.text(j,i,f'{piv.values[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white' if piv.values[i,j]>piv.values.max()*.6 else '#333')
 cb=fig.colorbar(im,ax=ax,fraction=0.04,pad=0.03); cb.ax.tick_params(labelsize=6)
save(fig,DEST/'S2_A_ECM_Notch_heatmaps')

# B: key ligand-receptor lollipop
fig,ax=plt.subplots(figsize=(5.8,4.4)); top=lr[lr.p_perm<0.05].copy().sort_values('prob',ascending=False).drop_duplicates(['pair','source','target']).head(8); top['label']=top.pair+'\n'+top.source+'→'+top.target; top=top.sort_values(['source','prob'],ascending=[True,True]); y=np.arange(len(top)); cols=['#E45756' if t=='PH' else '#72B7B2' for t in top.target]; ax.hlines(y,0,top.prob,color='#DDD',lw=1); ax.scatter(top.prob,y,c=cols,s=34,zorder=3)
for i,r in enumerate(top.itertuples()):
 star='***' if r.p_perm<.001 else '**' if r.p_perm<.01 else '*'; ax.text(r.prob+.004,i,star,va='center',fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(top.label,fontsize=6.3); ax.set_xlabel('Communication probability'); ax.set_xlim(0,top.prob.max()*1.18); ax.set_title('Key ligand–receptor pairs'); ax.legend(handles=[Patch(color='#E45756',label='PH'),Patch(color='#72B7B2',label='GH')],frameon=False,fontsize=6.5,loc='lower right')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,DEST/'S2_B_ligand_receptor_lollipop')

# C: ITGB1 violin
fig,ax=plt.subplots(figsize=(5.2,3.8)); groups={'PH':['PH'],'GH':['GH'],'Niche':['TS','MH2','MyoEpi','M'],'Hepato':['F','R','B'],'Other':['PM-like']}; data=[itgb1[a.obs.cell_type.isin(v).values] for v in groups.values()]; vp=ax.violinplot(data,showextrema=False,widths=.8)
for i,b in enumerate(vp['bodies']): b.set_facecolor('#C44E52' if i<2 else '#8DA0CB'); b.set_alpha(.72); b.set_edgecolor('#444'); b.set_linewidth(.4)
for i,d in enumerate(data): ax.plot([i,i],[np.median(d),np.median(d)],color='#222',lw=1.1)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(list(groups),rotation=30); ax.set_ylabel('ITGB1 (log1p CP10K)'); ax.set_title('ITGB1 receptor expression')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,DEST/'S2_C_ITGB1_violin')

# D: dynamic niche signalling volcano
fig,ax=plt.subplots(figsize=(5.5,4.2)); piv=gdf.groupby(['source','target','group'])['prob'].sum().unstack(fill_value=0).reset_index(); piv['label']=piv.source+'→'+piv.target; piv['log2fc']=np.log2((piv.Hemospasia+1e-6)/(piv.Control+1e-6)); labels=a.obs.cell_type.astype(str).values; grp=a.obs.group.astype(str).replace({'blank':'Control','treated':'Hemospasia'}).values; rng=np.random.default_rng(42); ps=[]
for r in piv.itertuples():
 pairs=lr[(lr.source==r.source)&(lr.target==r.target)]; lg=np.zeros(a.n_obs); rg=np.zeros(a.n_obs)
 for p_ in pairs.itertuples(): lg+=expr(str(p_.ligand_genes).split(';')); rg+=expr(str(p_.receptor_genes).split(';'))
 sm=labels==r.source; tm=labels==r.target; obs=lg[sm&(grp=='Hemospasia')].mean()*rg[tm&(grp=='Hemospasia')].mean()-lg[sm&(grp=='Control')].mean()*rg[tm&(grp=='Control')].mean(); null=[]
 for _ in range(500):
  perm=rng.permutation(grp); null.append(lg[sm&(perm=='Hemospasia')].mean()*rg[tm&(perm=='Hemospasia')].mean()-lg[sm&(perm=='Control')].mean()*rg[tm&(perm=='Control')].mean())
 null=np.array(null); ps.append((np.sum(np.abs(null)>=abs(obs))+1)/len(null))
piv['p']=ps; ax.scatter(piv.log2fc,-np.log10(piv.p.clip(lower=1e-4)),s=42,c=['#C44E52' if p<.1 else '#999' for p in piv.p],edgecolors='#333',linewidths=.3,zorder=3)
for r in piv.itertuples():
 if r.p<.1: ax.annotate(r.label,(r.log2fc,-np.log10(max(r.p,1e-4))),textcoords='offset points',xytext=(4,3),fontsize=6.2)
ax.axhline(-np.log10(.1),color='#999',ls='--',lw=.6); ax.axvline(0,color='#999',ls='--',lw=.6); ax.set_xlabel('log2 FC (Hemospasia / Control)'); ax.set_ylabel('−log10 P'); ax.set_title('Dynamic niche signalling')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,DEST/'S2_D_dynamic_niche_volcano')

# E: pseudotime ridge + compartment comparison
fig=plt.figure(figsize=(7.2,4.2)); gs=fig.add_gridspec(1,2,wspace=.38,left=.08,right=.98,bottom=.2,top=.85); obs=a.obs.copy(); obs['grp2']=obs.cell_type.astype(str).replace({'F':'Hepato','R':'Hepato','B':'Hepato','PM-like':'Other'}); cats=['PH','GH','MyoEpi','M','TS','MH2','Hepato','Other']; xs=np.linspace(0,.55,220)
ax=fig.add_subplot(gs[0,0])
for i,c in enumerate(cats[::-1]):
 v=obs.loc[obs.grp2==c,'projected_pseudotime'].dropna().values
 if len(v)<10: continue
 y=gaussian_kde(v)(xs); y=y/y.max()*.8; ax.fill_between(xs,i,i+y,color='#B2182B' if c=='PH' else ('#72B7B2' if c=='GH' else '#8DA0CB'),alpha=.75,lw=0); ax.text(.002,i+.12,c,ha='left',va='center',fontsize=6.4)
ax.set_yticks([]); ax.set_xlim(0,.5); ax.set_xlabel('Projected pseudotime'); ax.set_title('Pseudotime distribution by state')
for sp in ('top','right','left'): ax.spines[sp].set_visible(False)
ax=fig.add_subplot(gs[0,1]); data=[obs.loc[obs.grp2==c,'projected_pseudotime'].dropna().values for c in cats]; bp=ax.boxplot(data,tick_labels=cats,showfliers=False,patch_artist=True,widths=.6)
for p,c in zip(bp['boxes'],['#B2182B' if x=='PH' else ('#72B7B2' if x=='GH' else '#8DA0CB') for x in cats]): p.set_facecolor(c); p.set_alpha(.72)
ax.set_ylabel('Projected pseudotime'); ax.set_title('Compartment comparison'); ax.tick_params(axis='x',rotation=60)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig,DEST/'S2_E_pseudotime_ridge_box')
