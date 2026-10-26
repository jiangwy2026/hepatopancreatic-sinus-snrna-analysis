# -*- coding: utf-8 -*-
"""Export revised Figure 2A-C with larger axis fonts and corrected labels."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig1_Fig2_final" / "Figure2"
OUT.mkdir(parents=True, exist_ok=True)
CC = ROOT / "analysis_output" / "figures_pub" / "2026-10-21_Fig2" / "preview" / "Fig2a_all_states_cellcycle.csv"
COMP = ROOT / "analysis_output" / "Fig5_dynamic_subpop_subcluster_composition.csv"
HEAT = ROOT / "analysis_output" / "Fig5_PH_dpt_heatmap_data.csv"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family":"serif", "font.serif":["Times New Roman"], "font.size":8.5,
    "axes.titlesize":9.5, "axes.labelsize":9.5, "xtick.labelsize":8, "ytick.labelsize":8,
    "axes.linewidth":0.7, "svg.fonttype":"none", "pdf.fonttype":42, "savefig.facecolor":"white"
})

def save(fig, stem):
    fig.savefig(stem.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(stem.with_suffix('.svg'), format='svg', bbox_inches='tight')
    fig.savefig(stem.with_suffix('.pdf'), format='pdf', bbox_inches='tight')
    fig.savefig(stem.with_suffix('.tiff'), format='tiff', dpi=1200, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

# A: cell-cycle composition
d = pd.read_csv(CC).sort_values('G2M', ascending=True)
colors = {'G1':'#BDBDBD', 'S':'#72B7B2', 'G2M':'#B2182B'}
fig, ax = plt.subplots(figsize=(4.8,3.5))
left = np.zeros(len(d))
for phase in ['G1','S','G2M']:
    ax.barh(d.cell_type, d[phase], left=left, color=colors[phase], height=0.62, label=phase)
    left += d[phase].values
for i,r in enumerate(d.itertuples()):
    ax.text(101,i,f'G2M {r.G2M:.1f}% (n={r.n})',va='center',fontsize=6.8,fontweight='bold' if r.cell_type=='PH' else 'normal')
ax.set_xlim(0,145); ax.set_xlabel('Cells per phase (%)'); ax.set_ylabel('Cell state')
ax.set_title('Proliferation potential across ten Hs states', pad=5)
ax.legend(frameon=False,fontsize=7,ncol=3,loc='lower center',bbox_to_anchor=(0.5,-0.30))
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig, OUT/'Fig2A_proliferation')

# B: PH substate composition
d = pd.read_csv(COMP); d=d[d.population=='PH'].sort_values('subcluster')
labels={0:'s0 ECM/ribosome-high',1:'s1 progenitor-like',2:'s2 immune/ECM-primed',3:'s3 TGM1-low'}
cols={0:'#4C78A8',1:'#F58518',2:'#E45756',3:'#54A24B'}
fig,ax=plt.subplots(figsize=(4.5,3.3)); xpos=np.array([0.,0.78]); bottom=np.zeros(2)
for _,r in d.iterrows():
    sid=int(r.subcluster); vals=np.array([r.pct_control,r.pct_hemospasia])
    ax.bar(xpos,vals,bottom=bottom,color=cols[sid],width=0.60,edgecolor='white',linewidth=0.8,label=labels[sid])
    for i,(v,b) in enumerate(zip(vals,bottom)):
        if v>=6: ax.text(xpos[i],b+v/2,f's{sid} {v:.1f}%',ha='center',va='center',fontsize=6.8,color='white' if sid!=3 else '#222')
    bottom+=vals
ax.set_xticks(xpos); ax.set_xticklabels(['Control','Hemospasia']); ax.set_ylim(0,100); ax.set_ylabel('PH cells (%)')
ax.set_title('PH substate composition',pad=6)
ax.text(xpos[0],-8,f'n = {int(d.n_control.sum())}',ha='center',va='top',fontsize=7)
ax.text(xpos[1],-8,f'n = {int(d.n_hemospasia.sum())}',ha='center',va='top',fontsize=7)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.legend(handles=[Patch(color=cols[k],label=labels[k]) for k in sorted(labels)],frameon=False,fontsize=6.4,loc='center left',bbox_to_anchor=(1.02,0.5))
fig.tight_layout(); save(fig,OUT/'Fig2B_substate_composition')

# C: pseudotime-ordered marker heatmap
d=pd.read_csv(HEAT); genes=d.gene.drop_duplicates().tolist(); bins=sorted(d.bin.unique())
Z=np.zeros((len(genes),len(bins))); dom=[]
for i,g in enumerate(genes):
    sub=d[d.gene==g].set_index('bin').reindex(bins)
    Z[i]=sub.zscore.values
    dom.append(sub.dominant_substate.iloc[0] if len(sub) else 'PH_s1')
fig, ax = plt.subplots(figsize=(6.6,4.6))
fig.subplots_adjust(left=0.22, right=0.86, top=0.88, bottom=0.16)
im=ax.imshow(np.clip(Z,-2,2),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes,fontsize=7.8)
ax.set_xticks([0,len(bins)//2,len(bins)-1]); ax.set_xticklabels(['early','mid','late'])
ax.set_xlabel('Pseudotime',fontsize=9.5); ax.set_ylabel('PH marker genes',fontsize=9.5)
ax.set_title('PH marker dynamics along pseudotime',pad=8)
for sp in ax.spines.values(): sp.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.04,pad=0.03); cb.set_label('Scaled expression',fontsize=8); cb.ax.tick_params(labelsize=7)
save(fig,OUT/'Fig2C_pseudotime_heatmap')
print('Figure2 A-C saved')
