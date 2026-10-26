# -*- coding: utf-8 -*-
"""196_syncytium_summary_preview.py - summarize existing syncytium checks."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
os.makedirs(PREV,exist_ok=True)

def main():
    vals=pd.DataFrame({
        "set":["Retained cells","Retained + dropped PH-like","Dropped PH-like"],
        "blank":[6.84,4.19,1.97],"treated":[15.73,10.99,5.52],
        "fold":[2.30,2.62,2.80]})
    fig=plt.figure(figsize=(7.2,3.1),facecolor="white")
    ax=fig.add_axes([0.23,0.22,0.34,0.62]); y=np.arange(len(vals))
    ax.barh(y-0.17,vals.blank,height=0.32,color="#4C78A8",label="Control")
    ax.barh(y+0.17,vals.treated,height=0.32,color="#D6604D",label="Hemospasia")
    ax.set_yticks(y); ax.set_yticklabels(vals.set,fontsize=5.6); ax.invert_yaxis()
    ax.set_xlabel("PH proportion (%)",fontsize=5.8); ax.set_xlim(0,18)
    for i,r in vals.iterrows():
        ax.text(max(r.blank,r.treated)+0.5,i,"%.2f×"%r.fold,va="center",fontsize=5.0,color="#333333")
    ax.legend(frameon=False,fontsize=5.0,loc="lower right"); ax.set_title("PH expansion sensitivity",fontsize=6.8,pad=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x",labelsize=5.0,length=1.1); ax.tick_params(axis="y",length=0)

    # retained PH vs dropped PH-like signature
    markers=["PH signature","TGM","SSPO","ccnd2","Grik2","PPAF2","G2M score","pct_mt"]
    retained=[2.48,1.72,2.85,1.52,1.88,2.69,0.16,6.2]
    dropped=[2.60,1.86,2.93,1.61,2.10,2.88,0.19,15.08]
    Z=np.vstack([retained,dropped]).T
    # row normalize for display
    Z=Z/(np.max(Z,axis=1,keepdims=True)+1e-9)
    ax2=fig.add_axes([0.66,0.22,0.17,0.62])
    im=ax2.imshow(Z,aspect="auto",cmap="viridis",vmin=0,vmax=1)
    ax2.set_xticks([0,1]); ax2.set_xticklabels(["Retained\nPH","Dropped\nPH-like"],fontsize=5.0)
    ax2.set_yticks(range(len(markers))); ax2.set_yticklabels(markers,fontsize=5.0)
    ax2.tick_params(length=0)
    for i in range(len(markers)):
        for j in range(2):
            v=(retained if j==0 else dropped)[i]
            ax2.text(j,i,"%.2f"%v,ha="center",va="center",fontsize=4.4,color="white" if Z[i,j]<0.6 else "black")
    ax2.set_title("PH-like program",fontsize=6.2,pad=4)
    cb=fig.colorbar(im,ax=ax2,fraction=0.045,pad=0.03); cb.set_label("row-scaled",fontsize=4.8); cb.ax.tick_params(labelsize=4.4,length=0.7); cb.outline.set_visible(False)
    fig.text(0.50,0.035,"No intact 2× PH syncytium was detected; giant barcodes were mixed-cell clumps.",ha="center",va="center",fontsize=5.0,color="#333333")
    p=os.path.join(PREV,"07_syncytium_summary_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": main()
