# -*- coding: utf-8 -*-
"""227_lyo_cluster_mapping_bubble.py
LyO cluster marker signature AUROC bubble plot in the same style as the
previous hemocyte/hepatopancreas annotation bubbles.
"""
import os, sys, numpy as np, pandas as pd, scanpy as sc
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import rankdata
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','淋巴器官亚群映射')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42}); sc.settings.verbosity=0
NAMES={0:'ECM-producing',1:'Macrophage-like 1',2:'Macrophage-like 2',3:'Unknown 3',4:'Lymphocyte precursor-like',5:'Capsular',6:'Unknown 6',7:'Hemocytes',8:'SLC-rich',9:'Unknown 9',10:'Unknown 10',11:'Epithelial',12:'Unknown 12',13:'Progenitor',14:'Unknown 14',15:'Unknown 15'}
ROWS=['PH','GH','TS','MH2','F','R','B','M','MyoEpi','PM-like']
CIDS={'TS':'C0','GH':'C1','PH':'C2','MH2':'C3','PM-like':'C4','F':'C5','R':'C6','B':'C7','MyoEpi':'C8','M':'C9'}
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
def load_lyo():
 d=pd.read_excel(os.path.join(ROOT,'DeGene.list.xls'),sheet_name='Top50.DeGene.list')
 d['Target_Cluster']=d.Target_Cluster.astype(int)
 sig={}
 for c,g in d.groupby('Target_Cluster'):
  g=g.sort_values('Log2FC',ascending=False)
  sig[c]=list(dict.fromkeys(g['Gene Name'].astype(str).tolist()))[:50]
 return sig
def main():
 A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad')); A.obs['ct']=A.obs.vendor_cluster.astype(str).map(GROUPS)
 A.X=A.layers['counts'].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 sym2i={}
 for i,s in enumerate(A.var['sym'].astype(str)): sym2i.setdefault(s,i)
 sig=load_lyo(); rows=[]; ctvals=A.obs.ct.astype(str).values; n=A.n_obs
 for c in range(16):
  genes=[g for g in sig.get(c,[]) if g in sym2i]
  if len(genes)<5: continue
  X=np.asarray(A[:, [sym2i[g] for g in genes]].X.toarray(),dtype=float)
  pct=np.array([100*(X[ctvals==t]>0).mean(0) for t in ROWS])
  det=np.where((pct>1).any(0))[0]
  if len(det)==0: continue
  score=X[:,det].mean(1); r=rankdata(score)
  for i,ct in enumerate(ROWS):
   m=(ctvals==ct); npos=m.sum(); nneg=n-npos
   auc=(r[m].sum()-npos*(npos+1)/2)/(npos*nneg)
   other=np.delete(pct[:,det],i,axis=0); sup=int(np.sum((pct[i,det]>other.max(0)) & (pct[i,det]>5)))
   rows.append({'lyo_cluster':c,'lyo_name':NAMES[c],'target_group':ct,'n_genes':len(genes),'support_genes':sup,'auroc':float(auc)})
 R=pd.DataFrame(rows); R.to_csv(os.path.join(OUT,'Fig_LyO_mapping_AUROC.csv'),index=False,encoding='utf-8-sig')
 Z=R.pivot(index='target_group',columns='lyo_cluster',values='auroc').reindex(index=ROWS,columns=list(range(16)))
 N=R.pivot(index='target_group',columns='lyo_cluster',values='support_genes').reindex(index=ROWS,columns=list(range(16)))
 fig=plt.figure(figsize=(7.2,4.3),facecolor='white'); ax=fig.add_axes([0.11,0.27,0.68,0.60])
 cmap=plt.get_cmap('viridis'); norm=plt.Normalize(.5,1.0)
 for i,r in enumerate(ROWS):
  for j,c in enumerate(range(16)):
   a=Z.loc[r,c] if c in Z.columns else np.nan; n=N.loc[r,c] if c in N.columns else np.nan
   if pd.isna(a): ax.scatter(j,i,s=7,c='white',edgecolors='#CCCCCC',linewidths=.4); continue
   ax.scatter(j,i,s=10+3.4*n,c=[cmap(norm(a))],edgecolors='white',linewidths=.45,zorder=3)
 ax.set_xticks(range(16)); ax.set_xticklabels([NAMES[c] for c in range(16)],rotation=45,ha='right',rotation_mode='anchor',fontsize=4.7)
 ax.set_yticks(range(len(ROWS))); ax.set_yticklabels([f'{r}({CIDS[r]})' for r in ROWS],fontsize=5.2)
 ax.set_xlim(-.55,15.55); ax.set_ylim(len(ROWS)-.5,-.5); ax.tick_params(length=0,pad=2)
 for sp in ax.spines.values(): sp.set_visible(False)
 ax.grid(color='#F0F0F0',lw=.35); ax.set_axisbelow(True)
 cax=fig.add_axes([0.815,0.32,0.016,0.44]); cb=fig.colorbar(plt.cm.ScalarMappable(cmap=cmap,norm=norm),cax=cax); cb.set_label('AUROC',fontsize=5.2,labelpad=4); cb.ax.tick_params(labelsize=5.0,length=1.2,width=.4); cb.outline.set_linewidth(.4)
 axl=fig.add_axes([0.805,0.15,0.12,0.13]); axl.set_xlim(0,1); axl.set_ylim(0,1); axl.axis('off'); axl.text(0,.95,'genes',fontsize=5.0)
 for k,n in enumerate([5,10,20]): axl.scatter(.18,.68-k*.25,s=10+3.4*n,c='white',edgecolors='#666666',linewidths=.45); axl.text(.42,.68-k*.25,str(n),fontsize=5.0,va='center')
 fig.text(.5,.965,'Lymphoid-organ cluster signatures across hepatosinus cell states',ha='center',va='top',fontsize=7.2)
 fig.text(.5,.035,'Bubble size = supported LyO marker genes; colour = AUROC of the LyO signature score for the indicated cell state vs all other cells. Marker-level mapping, not label transfer.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig_LyO_mapping_bubble_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig_LyO_mapping_bubble.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig_LyO_mapping_bubble.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
