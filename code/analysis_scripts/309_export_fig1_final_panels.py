# -*- coding: utf-8 -*-
"""Export revised Figure 1 panels with unified tSNE axes and larger fonts."""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig1_Fig2_final" / "Figure1"
OUT.mkdir(parents=True, exist_ok=True)
TSNE_PATH = ROOT / "analysis_output" / "adata_soupx_hvg_tsne.h5ad"
A_PATH = ROOT / "analysis_output" / "adata_soupx.h5ad"
COMP_PATH = ROOT / "analysis_output" / "composition_new_names_10groups.csv"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})
GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
ORDER = ["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]
FULL = {"PH": "Prohemocyte", "GH": "Granular hemocyte", "TS": "Biosynthetic transition",
        "MH2": "Macrophage-like", "F": "F cell", "R": "R cell", "B": "B cell",
        "M": "M cell", "MyoEpi": "Myoepithelial-like", "PM-like": "Peritrophic matrix-like"}
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "M": "#2166AC",
       "MyoEpi": "#92C5DE", "PM-like": "#7B3294"}

T = sc.read_h5ad(TSNE_PATH)
A = sc.read_h5ad(A_PATH)[T.obs_names].copy()
ct = A.obs["vendor_cluster"].astype(str).map(GROUPS).values
sample = A.obs["sample"].astype(str).values
tsne = np.asarray(T.obsm["X_tsne"])
xmin, xmax = tsne[:,0].min(), tsne[:,0].max(); ymin, ymax = tsne[:,1].min(), tsne[:,1].max()
px, py = (xmax-xmin)*0.035, (ymax-ymin)*0.035
XLIM=(xmin-px,xmax+px); YLIM=(ymin-py,ymax+py)

def save(fig, stem):
    fig.savefig(stem.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(stem.with_suffix('.svg'), format='svg', bbox_inches='tight')
    fig.savefig(stem.with_suffix('.pdf'), format='pdf', bbox_inches='tight')
    fig.savefig(stem.with_suffix('.tiff'), format='tiff', dpi=1200, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

def style_tsne(ax):
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
    ax.set_xlabel('tSNE', fontsize=9.5); ax.set_ylabel('tSNE', fontsize=9.5)
    for sp in ax.spines.values():
        sp.set_color('#8F8F8F'); sp.set_linewidth(0.55)

def draw_tsne(ax, mask, label=True, point_size=1.6):
    for g in ORDER:
        m = mask & (ct == g)
        ax.scatter(tsne[m,0], tsne[m,1], s=point_size, c=COL[g], linewidths=0, rasterized=True)
    if label:
        for g in ORDER:
            m = mask & (ct == g)
            if m.sum() >= 20:
                x, y = np.median(tsne[m,0]), np.median(tsne[m,1])
                ax.text(x, y, g, ha='center', va='center', fontsize=6.5, weight='bold',
                        color='#1e293b', bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.78), zorder=8)

# B: QC
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
metrics = [('n_genes_by_counts','Genes per cell'), ('total_counts','UMIs per cell'), ('pct_counts_mt','Mitochondrial %')]
for ax,(col,label) in zip(axes, metrics):
    data=[]; pos=[]; cols=[]
    for j,grp in enumerate(['blank','treated']):
        v=A.obs.loc[sample==grp,col].astype(float).values
        data.append(v); pos.append(j); cols.append('#4C78A8' if grp=='blank' else '#E45756')
    bp=ax.violinplot(data, positions=pos, widths=0.65, showextrema=False)
    for body,c in zip(bp['bodies'],cols):
        body.set_facecolor(c); body.set_alpha(0.62); body.set_edgecolor('none')
    ax.scatter(pos,[np.median(x) for x in data],s=12,color='#222',zorder=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Control','Hemospasia'])
    ax.set_ylabel(label, fontsize=9.5); ax.set_title(label, fontsize=9)
    ax.tick_params(length=3,width=0.65)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig, OUT/'Fig1B_QC')

# C: all-cell tSNE
fig, ax = plt.subplots(figsize=(4.9,5.3))
fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.24)
draw_tsne(ax,np.ones(A.n_obs,dtype=bool),True,1.55)
style_tsne(ax); ax.set_title('All Hs cells (n=9,658)', pad=5)
handles=[Line2D([0],[0],marker='o',linestyle='',markersize=5,markerfacecolor=COL[g],markeredgewidth=0,label=f'{g} ({int((ct==g).sum()):,})') for g in ORDER]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5,0.015), ncol=2, frameon=False, fontsize=6.8, handletextpad=0.3, columnspacing=1.0)
save(fig,OUT/'Fig1C_tSNE_all')

# E/F: group split tSNE
for j,grp in enumerate(['blank','treated']):
    fig, ax = plt.subplots(figsize=(4.2,4.3))
    draw_tsne(ax,sample==grp,True,1.7)
    style_tsne(ax); ax.set_title(('Control (0 h)' if grp=='blank' else 'Hemospasia (6 h)')+f'  n={int((sample==grp).sum()):,}', pad=5)
    fig.tight_layout(); save(fig, OUT/f'Fig1{"E" if j==0 else "F"}_tSNE_{"control" if j==0 else "hemospasia"}')

# G: composition and fold change
d=pd.read_csv(COMP_PATH)
p=[]
for r in d.itertuples():
    _,pv=proportions_ztest([r.n_treated,r.n_blank],[4309,5349]); p.append(pv)
d['fdr']=multipletests(p,method='fdr_bh')[1]
d['log2fc']=np.log2(d.fold_treated_vs_blank.clip(lower=1e-6))
d=d.set_index('group').reindex(ORDER).reset_index()
fig=plt.figure(figsize=(7.4,3.2))
gs=fig.add_gridspec(1,2,width_ratios=[0.9,1.4],wspace=0.38,left=0.08,right=0.98,bottom=0.18,top=0.90)
ax=fig.add_subplot(gs[0,0])
x=np.arange(2); bottom=np.zeros(2)
for g in ORDER:
    r=d[d.group==g].iloc[0]
    vals=np.array([r.pct_blank,r.pct_treated])
    ax.bar(x,vals,bottom=bottom,width=0.62,color=COL[g],edgecolor='white',linewidth=0.5,label=g)
    bottom += vals
ax.set_xticks(x); ax.set_xticklabels(['Control','Hemospasia']); ax.set_ylim(0,100)
ax.set_ylabel('Cells (%)'); ax.set_title('Cell-state composition')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax=fig.add_subplot(gs[0,1])
ss=d.sort_values('log2fc')
y=np.arange(len(ss)); colors=['#C44E52' if v>0 else '#4C78A8' for v in ss.log2fc]
ax.hlines(y,0,ss.log2fc,color='#D0D0D0',lw=1.0)
ax.scatter(ss.log2fc,y,s=28,color=colors,zorder=3)
for i,r in enumerate(ss.itertuples()):
    ax.text(r.log2fc+(0.05 if r.log2fc>0 else -0.05),i,f'{r.fold_treated_vs_blank:.2f}×; FDR={r.fdr:.1e}',va='center',ha='left' if r.log2fc>0 else 'right',fontsize=6.2)
ax.axvline(0,color='#333',lw=0.6); ax.set_yticks(y); ax.set_yticklabels(ss.group)
ax.set_xlabel('log2 fold change (Hemospasia / Control)'); ax.set_title('Composition change')
ax.set_xlim(ss.log2fc.min()*1.7,ss.log2fc.max()*2.2)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.savefig(OUT/'Fig1G_composition.png',dpi=300,bbox_inches='tight'); fig.savefig(OUT/'Fig1G_composition.svg',format='svg',bbox_inches='tight'); fig.savefig(OUT/'Fig1G_composition.pdf',format='pdf',bbox_inches='tight'); fig.savefig(OUT/'Fig1G_composition.tiff',format='tiff',dpi=1200,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); plt.close(fig)
print('Figure1 panels saved')
