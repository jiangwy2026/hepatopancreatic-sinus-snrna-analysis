# -*- coding: utf-8 -*-
"""224_parent_program_effect_heatmap.py
QC-adjusted parent-state program effects (Hemospasia vs Control)."""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys, numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import t as tdist
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
sc.settings.verbosity=0
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=['PH','GH','TS','MH2','F','R','B','PM-like']
PROGS={'PH/progenitor':['PPAF2','SSPO','Grik2','ccnd2','TGM1','Syne1','slo','SERPINB1'],'Macrophage/immune':['TIE1','ITGA2B','PTPRA','SDC2','Tl','AGAP007347'],'Absorptive/metabolic':['ACE','Oacyl','HSPG2','mesh','hdc','NAAT1','SLC5A1','Npc1'],'ECM/remodeling':['THBS4','HSPG2','Tinagl1','Angpt4','SDC2'],'Stress':['hsc71','Hspa14','HSP90','HSPA','DnaJ'],'PM-like/cuticle':['Fas2','LOC113809538','LOC113809539','LOC113809540','GPx'],'G2M':['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']}
def bh(p):
 p=np.asarray(p,float); out=np.full_like(p,np.nan); ok=np.isfinite(p); idx=np.where(ok)[0]; order=idx[np.argsort(p[idx])]; m=len(order); prev=1.0
 for rank,i in enumerate(order[::-1],1):
  v=min(prev,p[i]*m/(m-rank+1)); out[i]=v; prev=v
 return out
def ols_group(score,group,nf,nc):
 y=np.asarray(score,float); g=np.asarray(group,float); nf=np.asarray(nf,float); nc=np.asarray(nc,float)
 nf=(nf-nf.mean())/(nf.std()+1e-9); nc=(nc-nc.mean())/(nc.std()+1e-9)
 X=np.column_stack([np.ones(len(y)),g,nf,nc]); beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
 resid=y-X@beta; df=len(y)-X.shape[1]; sd=np.sqrt(np.sum(resid**2)/df) if df>0 else np.nan
 inv=np.linalg.pinv(X.T@X); se=np.sqrt(max(sd**2*inv[1,1],1e-12)); tv=beta[1]/se; pv=2*tdist.sf(abs(tv),df) if df>0 else np.nan
 return beta[1], beta[1]/sd if sd>0 else np.nan, pv
def main():
 A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad')); A.X=A.layers['counts'].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 A.obs['parent']=A.obs.vendor_cluster.astype(str).map(GROUPS); A.obs['grp']=A.obs['sample'].astype(str).map({'blank':'Control','treated':'Hemospasia'}); A.obs['gind']=(A.obs.grp=='Hemospasia').astype(int)
 sym2gid={}
 for gid,s in zip(A.var_names.astype(str),A.var['sym'].astype(str)): sym2gid.setdefault(s,gid)
 for n,genes in PROGS.items(): sc.tl.score_genes(A,[sym2gid[g] for g in genes if g in sym2gid],score_name=n,random_state=0)
 rows=[]; pv=[]
 for par in ORDER:
  d=A.obs[A.obs.parent==par]
  for pr in PROGS:
   beta,dadj,p=ols_group(d[pr].values,d.gind.values,d.n_genes_by_counts.values,d.total_counts.values)
   c=d.loc[d.grp=='Control',pr]; t=d.loc[d.grp=='Hemospasia',pr]
   sd0=np.sqrt((np.var(c,ddof=1)+np.var(t,ddof=1))/2); d0=(t.mean()-c.mean())/sd0 if sd0>0 else 0
   rows.append({'parent':par,'program':pr,'n_control':len(c),'n_hemospasia':len(t),'mean_control':c.mean(),'mean_hemospasia':t.mean(),'cohens_d_unadjusted':d0,'cohens_d_qc_adjusted':dadj,'p_qc_adjusted':p}); pv.append(p)
 R=pd.DataFrame(rows); R['fdr_qc_adjusted']=bh(pv); R.to_csv(os.path.join(OUT,'Fig5_parent_program_effects_qc_adjusted.csv'),index=False,encoding='utf-8-sig')
 D=R.pivot(index='parent',columns='program',values='cohens_d_qc_adjusted').reindex(index=ORDER,columns=list(PROGS)); Q=R.pivot(index='parent',columns='program',values='fdr_qc_adjusted').reindex(index=ORDER,columns=list(PROGS))
 fig=plt.figure(figsize=(6.8,4.2),facecolor='white'); ax=fig.add_axes([0.18,0.26,0.68,0.62])
 im=ax.imshow(np.clip(D.values,-1.5,1.5),aspect='auto',cmap='RdBu_r',vmin=-1.5,vmax=1.5)
 ax.set_xticks(np.arange(len(PROGS))); ax.set_xticklabels(list(PROGS),rotation=38,ha='right',rotation_mode='anchor',fontsize=5.0)
 ax.set_yticks(np.arange(len(ORDER))); n={p:int((A.obs.parent==p).sum()) for p in ORDER}; ax.set_yticklabels([f"{p} (n={n[p]})" for p in ORDER],fontsize=5.2)
 ax.tick_params(length=0,pad=2)
 for i in range(len(ORDER)):
  for j in range(len(PROGS)):
   if Q.iloc[i,j]<0.05 and abs(D.iloc[i,j])>=0.2: ax.text(j,i,'*',ha='center',va='center',fontsize=6.5,color='#222222')
 for sp in ax.spines.values(): sp.set_visible(False)
 cax=fig.add_axes([0.875,0.30,0.018,0.44]); cb=fig.colorbar(im,cax=cax); cb.set_label("QC-adjusted Cohen's d",fontsize=5.2,labelpad=4); cb.ax.tick_params(labelsize=5.0,length=1.3,width=.4); cb.outline.set_linewidth(.4)
 fig.text(.5,.96,'QC-adjusted functional program changes',ha='center',va='top',fontsize=7.2)
 fig.text(.5,.03,"OLS adjusted for nFeature and nCount; * FDR < 0.05 and |d| >= 0.2. MyoEpi and M excluded by cell-number threshold.",ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_parent_program_effect_heatmap_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_parent_program_effects.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_parent_program_effects.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
