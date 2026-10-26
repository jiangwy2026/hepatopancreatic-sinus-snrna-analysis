# -*- coding: utf-8 -*-
"""286_ligand_specificity.py — ligand-side specificity for VWF/hemocytin, THBS4, POSTN."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig3_v2","preview"); os.makedirs(DEST,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,"axes.titlesize":8,
 "axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.linewidth":0.6,
 "savefig.facecolor":"white","svg.fonttype":"none"})
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad")); g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
lig={"VWF/hemocytin":["113824171"],"THBS4":["113823118"],"POSTN":["113817567","113805211"]}
rows=[]
for name,ids in lig.items():
    x=np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
    for ct in order:
        m=(a.obs.cell_type==ct).values
        rows.append({"ligand":name,"cell_type":ct,"mean":float(x[m].mean()),"pct":float((x[m]>0).mean()*100)})
d=pd.DataFrame(rows); d.to_csv(os.path.join(DEST,"Fig3_ligand_specificity.csv"),index=False,encoding="utf-8-sig")
print(d.pivot(index="ligand",columns="cell_type",values="pct").to_string())
print()
print(d.pivot(index="ligand",columns="cell_type",values="mean").to_string())
fig,ax=plt.subplots(figsize=(5.6,2.6))
for i,lig_name in enumerate(lig.keys()):
    for j,ct in enumerate(order):
        r=d[(d.ligand==lig_name)&(d.cell_type==ct)].iloc[0]
        ax.scatter(j,i,s=max(r.pct*1.6,2),c=r["mean"],cmap="Reds",vmin=0,vmax=d["mean"].max(),edgecolors="#666",linewidths=0.3)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=60)
ax.set_yticks(range(len(lig))); ax.set_yticklabels(list(lig))
for lab in ["MyoEpi","M"]:
    j=order.index(lab); ax.axvspan(j-0.5,j+0.5,color="#F4C7C3" if lab=="MyoEpi" else "#FCE3C3",alpha=0.25,zorder=0)
ax.set_title("Ligand-side specificity"); ax.set_xlim(-0.6,len(order)-0.4); ax.set_ylim(len(lig)-0.5,-0.5)
for v in [0,25,50,75]:
    ax.scatter([],[],s=max(v*1.6,2),c="#BBB",edgecolors="#666",linewidths=0.3,label=f"{v}%")
ax.legend(title="% expressing",frameon=False,fontsize=5.2,title_fontsize=5.2,loc="center left",bbox_to_anchor=(1.01,0.5))
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout()
for ext,kw in [("png",dict(dpi=300)),("svg",dict(format="svg")),("tiff",dict(format="tiff",dpi=1200,pil_kwargs={"compression":"tiff_lzw"}))]:
    fig.savefig(os.path.join(DEST,"Fig3_ligand_specificity."+ext),bbox_inches="tight",**kw)
print("saved")
