# -*- coding: utf-8 -*-
"""200_focused_lr_preview.py - top candidate ligand-receptor interactions."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
R=pd.read_csv(os.path.join(OUT,"focused_LR_candidates.csv"))
R["label"]=R.ligand+" ("+R.sender+") → "+R.receptor+" ("+R.receiver+")"
R=R.sort_values("score",ascending=False).head(20).iloc[::-1]
cols={"ECM-syndecan":"#8C510A","VEGF/PDGF":"#2166AC","FGF":"#41AB5D","Wnt":"#7B3294","IGF":"#D6604D","ECM-integrin":"#B2182B"}
fig=plt.figure(figsize=(5.6,4.0),facecolor="white")
ax=fig.add_axes([0.44,0.14,0.54,0.78])
ax.barh(np.arange(len(R)),R.score,color=[cols.get(x,"#777777") for x in R.pathway],height=0.65)
ax.set_yticks(np.arange(len(R))); ax.set_yticklabels(R.label,fontsize=4.8)
ax.set_xlabel("Candidate interaction score",fontsize=5.5); ax.set_title("Focused ligand–receptor screen",fontsize=6.8,pad=5)
ax.tick_params(axis="x",labelsize=4.8,length=1.0); ax.tick_params(axis="y",length=0)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.text(0.50,0.025,"Orthology-based candidate screen; not a validated shrimp ligand–receptor database.",ha="center",va="center",fontsize=4.7,color="#555555")
p=os.path.join(PREV,"10_focused_LR_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)
