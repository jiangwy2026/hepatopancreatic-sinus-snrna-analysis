# -*- coding: utf-8 -*-
"""Export Figure 2F: Hs PH versus hemolymph PH, with corrected naming."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

ROOT = Path(r"D:\桌面\single cell")
SRC = ROOT / "analysis_output" / "cross_organ_progenitor"
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig1_Fig2_final" / "Figure2"
OUT.mkdir(parents=True, exist_ok=True)
perc = pd.read_csv(SRC / "cross_organ_percell_scores.csv.gz")
strat = pd.read_csv(SRC / "cross_organ_cycle_stratified_means.csv")
summary = pd.read_csv(SRC / "cross_organ_adjusted_summary.csv")
order = ["Hs PH", "Blood PH"]
labels = {"Hs PH":"Hs PH", "Blood PH":"Hemolymph PH"}
colors = {"Hs PH":"#B2182B", "Blood PH":"#2166AC"}
perc = perc[perc.target.isin(order)].copy(); strat = strat[strat.target.isin(order)].copy(); summary = summary[summary.target.isin(order)].copy()

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family":"serif", "font.serif":["Times New Roman"], "font.size":8.5,
    "axes.titlesize":9.5, "axes.labelsize":9.5, "xtick.labelsize":7.8, "ytick.labelsize":7.8,
    "axes.linewidth":0.7, "svg.fonttype":"none", "pdf.fonttype":42, "savefig.facecolor":"white"
})

fig = plt.figure(figsize=(8.0,4.7))
gs = fig.add_gridspec(1,3,width_ratios=[1.0,1.0,0.95],wspace=0.42,left=0.07,right=0.98,bottom=0.17,top=0.88)
for j,key in enumerate(order):
    ax=fig.add_subplot(gs[0,j]); d=perc[perc.target==key]; bg=d[~d.is_target]; tg=d[d.is_target]
    ax.hexbin(bg.cycle_pct,bg.progenitor_pct,gridsize=32,cmap='Greys',mincnt=1,linewidths=0)
    ax.scatter(tg.cycle_pct,tg.progenitor_pct,s=8,alpha=0.28,color=colors[key],edgecolors='none',rasterized=True)
    sm=strat[(strat.target==key)&(strat.group=='background')]
    ax.plot(sm.mean_cycle_pct,sm.mean_prog_pct,color=colors[key],lw=1.3)
    row=summary[summary.target==key].iloc[0]
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xticks(np.linspace(0,1,5)); ax.set_yticks(np.linspace(0,1,5))
    ax.set_xlabel('Cell-cycle percentile'); ax.set_ylabel('Progenitor program percentile')
    ax.set_title(f"{labels[key]} (n={int(row.n_target):,})",pad=5)
    ax.tick_params(length=3,width=0.65)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

ax=fig.add_subplot(gs[0,2])
positions=np.arange(2)
for i,key in enumerate(order):
    d=perc[perc.target==key]; bgv=d.loc[~d.is_target,'adjusted_progenitor_z'].values; tgv=d.loc[d.is_target,'adjusted_progenitor_z'].values
    vp1=ax.violinplot([bgv],positions=[i-0.22],widths=0.42,bw_method=0.35,points=220,showextrema=False)
    vp2=ax.violinplot([tgv],positions=[i+0.22],widths=0.42,bw_method=0.35,points=220,showextrema=False)
    for body in vp1['bodies']: body.set_facecolor('#BDBDBD'); body.set_edgecolor('none'); body.set_alpha(0.72)
    for body in vp2['bodies']: body.set_facecolor(colors[key]); body.set_edgecolor('none'); body.set_alpha(0.78)
    ax.plot([i-0.22,i+0.22],[np.median(bgv),np.median(tgv)],color='#333',lw=0.8,zorder=3)
    ax.scatter([i-0.22,i+0.22],[np.median(bgv),np.median(tgv)],s=10,color=['#777',colors[key]],zorder=4)
    row=summary[summary.target==key].iloc[0]
    ax.text(i+0.22,0.82,f"z={row.mean_adjusted_prog_z:.2f}",transform=ax.get_xaxis_transform(),ha='center',va='top',color=colors[key],fontsize=7.2)
ax.axhline(0,color='#777',lw=0.7,ls='--'); ax.set_ylim(-2.6,2.6)
ax.set_xticks(positions); ax.set_xticklabels([labels[k] for k in order],fontsize=8)
ax.set_ylabel('Adjusted progenitor potential (z)'); ax.set_title('Progenitor identity after cell-cycle adjustment',pad=6)
ax.tick_params(length=3,width=0.65)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.legend(handles=[Patch(facecolor='#BDBDBD',label='background'),Patch(facecolor=colors['Hs PH'],label='Hs PH target'),Patch(facecolor=colors['Blood PH'],label='Hemolymph PH target')],frameon=False,fontsize=6.2,loc='upper left',ncol=1)

stem=OUT/'Fig2F_HsPH_vs_HemolymphPH'
fig.savefig(stem.with_suffix('.png'),dpi=300,bbox_inches='tight')
fig.savefig(stem.with_suffix('.svg'),format='svg',bbox_inches='tight')
fig.savefig(stem.with_suffix('.pdf'),format='pdf',bbox_inches='tight')
fig.savefig(stem.with_suffix('.tiff'),format='tiff',dpi=1200,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print(stem)
