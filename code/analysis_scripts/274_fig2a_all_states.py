# -*- coding: utf-8 -*-
"""274_fig2a_all_states.py — G2M potential across all 10 states."""
import os, collections
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"
DEST=os.path.join(OUT,"figures_pub","2026-10-21_Fig2","preview")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7.5,
 "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
 "axes.linewidth":0.6,"savefig.facecolor":"white","svg.fonttype":"none"})

a=ad.read_h5ad(os.path.join(OUT,"step3","step3_lyo_scores.h5ad"))
a.X=a.layers["lognorm"].copy()
g2v={str(s).lower(): a.var_names[i] for i,s in enumerate(a.var["sym"].astype(str))}
def pick(names):
    out=[]
    for n in names:
        v=g2v.get(n.lower())
        if v is not None and v not in out: out.append(v)
    return out
S=['PCNA','MCM6','CDC6','ORC1','ORC4','ORC5','POLA2','PRIM1','FEN1','TYMS','DHFR','HIRA','CHAF1A']
G2M=['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']
sg=pick(S); gg=pick(G2M)
print("S genes",len(sg),sg); print("G2M genes",len(gg),gg)
sc.tl.score_genes_cell_cycle(a,s_genes=sg,g2m_genes=gg)
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
rows=[]
for ct in order:
    m=(a.obs.cell_type==ct).values
    ph=a.obs.loc[m,"phase"].value_counts(normalize=True)*100
    rows.append({"cell_type":ct,"n":int(m.sum()),"G1":ph.get("G1",0),"S":ph.get("S",0),"G2M":ph.get("G2M",0)})
d=pd.DataFrame(rows).sort_values("G2M",ascending=True)
print(d.to_string(index=False))
d.to_csv(os.path.join(DEST,"Fig2a_all_states_cellcycle.csv"),index=False,encoding="utf-8-sig")

fig,ax=plt.subplots(figsize=(4.2,3.0))
left=np.zeros(len(d))
cols={"G1":"#BAB0AC","S":"#72B7B2","G2M":"#B2182B"}
for ph,col in cols.items():
    ax.barh(d.cell_type,d[ph],left=left,color=col,label=ph,height=0.62); left+=d[ph].values
for i,r in enumerate(d.itertuples()):
    ax.text(101,i,f"G2M {r.G2M:.1f}%  (n={r.n})",va="center",fontsize=5.8,fontweight="bold" if r.cell_type=="PH" else "normal")
ax.set_xlim(0,145);
for i,r in enumerate(d.itertuples()):
    if r.cell_type=="PH": ax.axhspan(i-0.4,i+0.4,color="#F4C7C3",alpha=0.35,zorder=0)
ax.set_xlabel("Cells per phase (%)")
ax.set_title("Proliferation potential across 10 states")
ax.legend(frameon=False,fontsize=6,ncol=3,loc="lower center",bbox_to_anchor=(0.5,-0.3))
ax.text(0.5,-0.42,"M: n=45 (small-population estimate); PH is the dominant proliferative state",transform=ax.transAxes,ha="center",fontsize=5.5,color="#666")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(DEST,"Fig2a_cellcycle_v3.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(DEST,"Fig2a_cellcycle_v3.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"Fig2a_cellcycle_v3.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
print("saved")
