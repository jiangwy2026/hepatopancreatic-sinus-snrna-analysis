# -*- coding: utf-8 -*-
"""229_pathway_enrichment_summary.py
Hypergeometric enrichment of DEGs against downloaded GO/KEGG/Reactome GMTs."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,'analysis_output'); PARTS=os.path.join(OUT,'figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化'); DE=os.path.join(OUT,'dynamic_subpop_de'); GS=os.path.join(OUT,'genesets')
os.makedirs(PARTS,exist_ok=True); os.makedirs(DELIV,exist_ok=True)
sys.path.insert(0,OUT)
from pub_style import apply_pub_style
apply_pub_style(7); mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
PARENTS=['PH','MH2','R','PM-like']; LIBS={'GO_BP':'GO_BP_2023.gmt','KEGG':'KEGG_2021_Human.gmt','Reactome':'Reactome_2022.gmt'}
def load_gmt(path):
 out=[]
 for line in open(path,encoding='utf-8',errors='replace'):
  parts=line.rstrip('\n').split('\t'); term=parts[0]; genes=set(x.upper() for x in parts[2:] if x)
  if 5<=len(genes)<=1000: out.append((term,genes))
 return out
def bh(p):
 p=np.asarray(p,float); out=np.full_like(p,np.nan); ok=np.isfinite(p); idx=np.where(ok)[0]; order=idx[np.argsort(p[idx])]; m=len(order); prev=1.0
 for rank,i in enumerate(order[::-1],1):
  v=min(prev,p[i]*m/(m-rank+1)); out[i]=v; prev=v
 return out
def main():
 import scanpy as sc
 A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad'),backed='r'); syms=set(str(x).upper() for x in A.var['sym'].astype(str)); A.file.close()
 bg={s for s in syms if not s.startswith(('LOC','GENE-'))}
 libs={k:load_gmt(os.path.join(GS,v)) for k,v in LIBS.items()}
 rows=[]
 for par in PARENTS:
  d=pd.read_csv(os.path.join(DE,f'DE_{par}_treated_vs_control.csv'))
  d=d[(d.pvals_adj<0.05)&(d.logfoldchanges.abs()>0.5)]
  for direction in ['up','down']:
   q=set(str(x).upper() for x in d.loc[(d.logfoldchanges>0) if direction=='up' else (d.logfoldchanges<0),'sym'].dropna())
   q={x for x in q if x in bg and not x.startswith(('LOC','GENE-'))}
   if len(q)<5: continue
   for libname,sets in libs.items():
    pv=[]; rec=[]
    for term,gs in sets:
     ov=len(q&gs)
     if ov<3: continue
     p=hypergeom.sf(ov-1,len(bg),len(gs),len(q)); pv.append(p); rec.append((term,gs,ov,p))
    if not pv: continue
    fd=bh(pv)
    for (term,gs,ov,p),f in zip(rec,fd):
     rows.append({'parent':par,'direction':direction,'library':libname,'term':term,'overlap':ov,'query_size':len(q),'term_size':len(gs),'p':p,'fdr':f,'genes':';'.join(sorted(q&gs))})
 R=pd.DataFrame(rows); R.to_csv(os.path.join(OUT,'Fig5_pathway_enrichment.csv'),index=False,encoding='utf-8-sig')
 cols=[f'{p}_{d}' for p in PARENTS for d in ['up','down']]; selected=[]; seen=set()
 for col in cols:
  par,di=col.rsplit('_',1); sub=R[(R.parent==par)&(R.direction==di)&(R.fdr<0.1)].sort_values(['fdr','overlap'],ascending=[True,False]).head(3)
  for r in sub.itertuples():
   key=f'{r.library}:{r.term}'
   if key not in seen: selected.append(key); seen.add(key)
 if not selected: print('No enriched terms'); return
 mat=np.full((len(selected),len(cols)),np.nan); ovm=np.zeros_like(mat)
 for i,key in enumerate(selected):
  lib,term=key.split(':',1)
  for j,col in enumerate(cols):
   par,di=col.rsplit('_',1); sub=R[(R.parent==par)&(R.direction==di)&(R.library==lib)&(R.term==term)]
   if len(sub): mat[i,j]=-np.log10(max(float(sub.iloc[0].fdr),1e-10)); ovm[i,j]=float(sub.iloc[0].overlap)
 fig=plt.figure(figsize=(7.2,6.4),facecolor='white'); ax=fig.add_axes([0.39,0.17,0.44,0.68])
 for i in range(len(selected)):
  for j in range(len(cols)):
   if np.isfinite(mat[i,j]):
    ax.scatter(j,i,s=10+5*ovm[i,j],c=[plt.cm.viridis(np.clip(mat[i,j]/8,0,1))],edgecolors='white',linewidths=.35,zorder=3)
 ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace('_',' ') for c in cols],rotation=45,ha='right',rotation_mode='anchor',fontsize=4.8)
 labels=[]
 for key in selected:
  lib,term=key.split(':',1); labels.append(('' if lib=='GO_BP' else lib+' ') + (term[:52]+'…' if len(term)>52 else term))
 ax.set_yticks(range(len(selected))); ax.set_yticklabels(labels,fontsize=4.7)
 ax.set_xlim(-.6,len(cols)-.4); ax.set_ylim(len(selected)-.5,-.5); ax.tick_params(length=0,pad=2)
 for sp in ax.spines.values(): sp.set_visible(False)
 ax.grid(color='#F2F2F2',lw=.35); ax.set_axisbelow(True)
 cax=fig.add_axes([0.86,0.30,0.017,0.42]); cb=fig.colorbar(plt.cm.ScalarMappable(cmap='viridis',norm=plt.Normalize(0,8)),cax=cax); cb.set_label('−log10 FDR',fontsize=5.2,labelpad=4); cb.ax.tick_params(labelsize=5.0,length=1.2,width=.4); cb.outline.set_linewidth(.4)
 axl=fig.add_axes([0.855,0.15,0.12,0.12]); axl.set_xlim(0,1); axl.set_ylim(0,1); axl.axis('off'); axl.text(0,.95,'Overlap',fontsize=5.0)
 for k,o in enumerate([5,10,20]): axl.scatter(.18,.68-k*.28,s=10+5*o,c='white',edgecolors='#666666',linewidths=.45); axl.text(.4,.68-k*.28,str(o),fontsize=5.0,va='center')
 fig.text(.5,.965,'Pathway enrichment of treatment-associated DEGs',ha='center',va='top',fontsize=7.2)
 fig.text(.5,.03,'Hypergeometric test against GO/KEGG/Reactome gene sets; LOC/unnamed genes excluded. Top 3 terms per column shown.',ha='center',va='bottom',fontsize=5.0,color='#555555')
 stem=os.path.join(PARTS,'Fig5_pathway_enrichment_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
 fig.savefig(os.path.join(DELIV,'Fig5_pathway_enrichment.svg'),bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(os.path.join(DELIV,'Fig5_pathway_enrichment.tiff'),dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 print('saved',stem); plt.close(fig)
if __name__=='__main__': main()
