# -*- coding: utf-8 -*-
"""223_qc_standardized_composition.py
Direct QC standardization of cell-state composition by nFeature/nCount strata.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-20", "第二张图_抽血处理细胞变化")
os.makedirs(PARTS, exist_ok=True); os.makedirs(DELIV, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","Nimbus Roman","Liberation Serif","DejaVu Serif"],"svg.fonttype":"none","pdf.fonttype":42})
sc.settings.verbosity=0

GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}

def main():
    A=sc.read_h5ad(os.path.join(OUT,'adata_soupx.h5ad'))
    A.obs['ct']=A.obs['vendor_cluster'].astype(str).map(GROUPS)
    A.obs['grp']=A.obs['sample'].astype(str).map({'blank':'Control','treated':'Hemospasia'})
    df=A.obs[['ct','grp','n_genes_by_counts','total_counts']].copy()
    # pooled QC strata
    df['nb']=pd.qcut(df.n_genes_by_counts,5,labels=False,duplicates='drop')
    df['cb']=pd.qcut(df.total_counts,5,labels=False,duplicates='drop')
    df['stratum']=df.nb.astype(str)+'_'+df.cb.astype(str)
    # pooled stratum weights
    w=df.stratum.value_counts(normalize=True).to_dict()
    groups=['Control','Hemospasia']; cats=sorted(df.ct.unique())
    raw={}; adj={}
    for g in groups:
        d=df[df.grp==g]
        raw[g]=d.ct.value_counts(normalize=True).reindex(cats,fill_value=0)
        vals={c:0.0 for c in cats}; wsum=0.0
        for s,ws in w.items():
            ds=d[d.stratum==s]
            if len(ds)==0: continue
            p=ds.ct.value_counts(normalize=True).reindex(cats,fill_value=0)
            for c in cats: vals[c]+=ws*p[c]
            wsum+=ws
        adj[g]=pd.Series({c: vals[c]/wsum if wsum>0 else np.nan for c in cats})
    rows=[]
    for c in cats:
        rp_c=100*raw['Control'][c]; rp_t=100*raw['Hemospasia'][c]
        ap_c=100*adj['Control'][c]; ap_t=100*adj['Hemospasia'][c]
        rows.append({'cluster':c,'n_control':int((df[(df.grp=='Control')&(df.ct==c)]).shape[0]),'n_hemospasia':int((df[(df.grp=='Hemospasia')&(df.ct==c)]).shape[0]),
                     'pct_control_raw':rp_c,'pct_hemospasia_raw':rp_t,
                     'fold_raw':rp_t/rp_c if rp_c>0 else np.nan,
                     'log2fc_raw':np.log2(rp_t/rp_c) if rp_c>0 and rp_t>0 else np.nan,
                     'pct_control_qc_std':ap_c,'pct_hemospasia_qc_std':ap_t,
                     'fold_qc_std':ap_t/ap_c if ap_c>0 else np.nan,
                     'log2fc_qc_std':np.log2(ap_t/ap_c) if ap_c>0 and ap_t>0 else np.nan})
    R=pd.DataFrame(rows).sort_values('log2fc_raw')
    R.to_csv(os.path.join(OUT,'Fig5_qc_standardized_composition.csv'),index=False,encoding='utf-8-sig')
    fig=plt.figure(figsize=(4.8,3.7),facecolor='white')
    ax=fig.add_axes([0.24,0.18,0.66,0.72])
    y=np.arange(len(R))
    for i,r in enumerate(R.itertuples()):
        ax.plot([r.log2fc_raw,r.log2fc_qc_std],[i,i],color='#BBBBBB',lw=1.0,zorder=1)
    ax.scatter(R.log2fc_raw,y,s=16,c='#777777',edgecolors='white',linewidths=.3,zorder=3,label='Unadjusted')
    ax.scatter(R.log2fc_qc_std,y,s=18,c='#B2182B',edgecolors='white',linewidths=.3,zorder=4,label='QC-standardized')
    ax.axvline(0,color='#333333',lw=.7)
    ax.set_yticks(y); ax.set_yticklabels([f"{c} ({int(r.n_control)}/{int(r.n_hemospasia)})" for c,r in zip(R.cluster,R.itertuples())],fontsize=5.2)
    ax.set_xlabel('log$_2$ fold change (Hemospasia / Control)',fontsize=5.6)
    ax.tick_params(axis='x',labelsize=5.0,length=1.4,width=.4); ax.tick_params(axis='y',length=0,pad=2)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    ax.grid(axis='x',color='#EEEEEE',lw=.35,zorder=0)
    ax.legend(loc='lower right',frameon=False,fontsize=5.0,handlelength=.9,handletextpad=.3)
    ax.set_title('Composition before and after QC standardization',fontsize=6.8,pad=5)
    fig.text(.5,.025,'QC strata: quintiles of nFeature and nCount; standardized to the pooled QC distribution. Cell-level descriptive analysis.',ha='center',va='bottom',fontsize=5.0,color='#555555')
    stem=os.path.join(PARTS,'Fig5_qc_standardized_composition_preview.png'); fig.savefig(stem,dpi=450,bbox_inches='tight',pad_inches=.04,facecolor='white')
    svg=os.path.join(DELIV,'Fig5_composition_QC_standardized.svg'); tif=os.path.join(DELIV,'Fig5_composition_QC_standardized.tiff')
    fig.savefig(svg,bbox_inches='tight',pad_inches=.04,facecolor='white'); fig.savefig(tif,dpi=1200,bbox_inches='tight',pad_inches=.04,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
    print('saved',stem); print('saved',svg); print('saved',tif)
    plt.close(fig)

if __name__=='__main__': main()
