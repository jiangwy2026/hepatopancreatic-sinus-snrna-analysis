# -*- coding: utf-8 -*-
"""230_coexpression_modules.py
Simple WGCNA-like co-expression modules and module-group associations."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import numpy as np, pandas as pd, scanpy as sc
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
sc.settings.verbosity=0
MODULE_LABELS={"M1":"PM-like/cuticle","M2":"Apod/ECM-lipid","M3":"absorptive/metabolic","M4":"mixed/uncertain","M5":"ECM/contractile","M6":"transport/membrane","M7":"PH progenitor"}
G={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}; PARENTS=['PH','GH','TS','MH2','F','R','B','PM-like']
def bh(p):
 p=np.asarray(p,float); out=np.full_like(p,np.nan); ok=np.isfinite(p); idx=np.where(ok)[0]; order=idx[np.argsort(p[idx])]; m=len(order); prev=1.0
 for rank,i in enumerate(order[::-1],1):
  v=min(prev,p[i]*m/(m-rank+1)); out[i]=v; prev=v
 return out
def main():
 A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad')); A.X=A.layers['counts'].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 sc.pp.highly_variable_genes(A,n_top_genes=1200,flavor='seurat'); H=A[:,A.var.highly_variable].copy(); X=np.asarray(H.X.toarray(),dtype=np.float32)
 C=np.corrcoef(X,rowvar=False); C=np.nan_to_num(C,nan=0,posinf=0,neginf=0); np.fill_diagonal(C,1)
 Z=linkage(squareform(1-C,checks=False),method='average'); mod=fcluster(Z,10,criterion='maxclust')
 genes=H.var['sym'].astype(str).values; M=pd.DataFrame({'gene':genes,'gid':H.var_names.astype(str),'module':mod})
 M.to_csv(os.path.join(OUT,'Fig5_coexpression_modules.csv'),index=False,encoding='utf-8-sig')
 scores=pd.DataFrame(index=A.obs_names)
 for k in sorted(set(mod)):
  idx=np.where(mod==k)[0]
  if len(idx)<5: continue
  Xs=X[:,idx]; Xs=(Xs-Xs.mean(0))/(Xs.std(0)+1e-6); pc=PCA(n_components=1,random_state=0).fit_transform(Xs).ravel()
  if np.corrcoef(pc,Xs.mean(1))[0,1]<0: pc=-pc
  scores[f'M{k}']=pc
 A.obs=A.obs.join(scores)
 A.obs['parent']=A.obs.vendor_cluster.astype(str).map(G); A.obs['grp']=A.obs['sample'].astype(str).map({'blank':'Control','treated':'Hemospasia'})
 rows=[]; pv=[]
 for par in PARENTS:
  d=A.obs[A.obs.parent==par]
  for mc in scores.columns:
   c=d.loc[d.grp=='Control',mc].values; t=d.loc[d.grp=='Hemospasia',mc].values
   sd=np.sqrt((np.var(c,ddof=1)+np.var(t,ddof=1))/2); dd=(np.mean(t)-np.mean(c))/sd if sd>0 else 0
   p=mannwhitneyu(c,t,alternative='two-sided').pvalue
   rows.append({'parent':par,'module':mc,'n_control':len(c),'n_hemospasia':len(t),'cohens_d':dd,'p':p}); pv.append(p)
 R=pd.DataFrame(rows); R['fdr']=bh(pv); R.to_csv(os.path.join(OUT,'Fig5_coexpression_module_effects.csv'),index=False,encoding='utf-8-sig')
 D=R.pivot(index='module',columns='parent',values='cohens_d').reindex(index=scores.columns,columns=PARENTS); Q=R.pivot(index='module',columns='parent',values='fdr').reindex(index=scores.columns,columns=PARENTS)
 fig=plt.figure(figsize=(6.8,4.2),facecolor='white'); ax=fig.add_axes([0.16,0.22,0.70,0.66])
 im=ax.imshow(np.clip(D.values,-1.5,1.5),aspect='auto',cmap='RdBu_r',vmin=-1.5,vmax=1.5)
 ax.set_xticks(range(len(PARENTS))); ax.set_xticklabels(PARENTS,rotation=38,ha='right',rotation_mode='anchor',fontsize=5.0)
 ax.set_yticks(range(len(scores.columns))); ax.set_yticklabels([f'{mc} {MODULE_LABELS.get(mc, mc)}' for mc in scores.columns],fontsize=5.2); ax.tick_params(length=0,pad=2)
 for i in range(len(scores.columns)):
  for j in range(len(PARENTS)):
   if Q.iloc[i,j]<0.05 and abs(D.iloc[i,j])>=0.2: ax.text(j,i,'*',ha='center',va='center',fontsize=6)
 for sp in ax.spines.values(): sp.set_visible(False)
 cax=fig.add_axes([0.875,0.30,0.018,0.42]); cb=fig.colorbar(im,cax=cax); cb.set_label("Cohen's d",fontsize=5.2,labelpad=4); cb.ax.tick_params(labelsize=5.0,length=1.2,width=.4); cb.outline.set_linewidth(.4)
 fig.text(.5,.96,'Co-expression module activity changes',ha='center',va='top',fontsize=7.2)
 fig.text(.5,.03,'Modules from 1,200 HVGs (WGCNA-like average-linkage clustering). * FDR < 0.05 and |d| >= 0.2.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_coexpression_module_heatmap_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_coexpression_module_heatmap.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_coexpression_module_heatmap.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
