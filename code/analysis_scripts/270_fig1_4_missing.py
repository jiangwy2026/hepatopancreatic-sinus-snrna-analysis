# -*- coding: utf-8 -*-
"""270_fig1_4_missing.py — missing panels: Fig1 composition stats, Fig4 ligand group, Fig4 hepatopancreas composition."""
import os
import numpy as np, pandas as pd, anndata as ad
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Figure_missing")
for s in ["svg","tiff","png"]: os.makedirs(os.path.join(DEST,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def save(fig,n):
    fig.savefig(os.path.join(DEST,"svg",n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",n+".png"),dpi=220,bbox_inches="tight"); plt.close(fig); print("saved",n,flush=True)

# Fig1 composition stats
d=pd.read_csv(os.path.join(OUT,"composition_new_names_10groups.csv"))
d["log2fc"]=np.log2(d["fold_treated_vs_blank"].clip(lower=1e-6))
p=[]
for r in d.itertuples():
    _,pv=proportions_ztest([r.n_treated,r.n_blank],[4309,5349]); p.append(pv)
d["p"]=p; d["fdr"]=multipletests(p,method="fdr_bh")[1]
d.to_csv(os.path.join(DEST,"Fig1_composition_stats.csv"),index=False,encoding="utf-8-sig")
fig,ax=plt.subplots(figsize=(4.2,3.0))
d=d.sort_values("log2fc")
ax.barh(d.group,d.log2fc,color=["#C44E52" if x>0 else "#4C78A8" for x in d.log2fc],height=0.62)
for i,r in enumerate(d.itertuples()):
    ax.text(r.log2fc+(0.03 if r.log2fc>0 else -0.03),i,f"{r.fold_treated_vs_blank:.2f}×; FDR={r.fdr:.1e}",
            va="center",ha="left" if r.log2fc>0 else "right",fontsize=4.8)
ax.axvline(0,color="#333",lw=0.6); ax.set_xlabel("log2 fold change (Hemospasia / Control)")
ax.set_title("Composition change (n=5,349 vs 4,309)"); ax.set_xlim(d.log2fc.min()*1.6,d.log2fc.max()*2.0)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); save(fig,"Fig1g_composition_stats")

# Fig4 ligand by group
a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids):
    idx=[g2v[i] for i in ids if i in g2v]; return np.asarray(a[:,idx].layers["lognorm"].todense()).mean(1) if idx else np.zeros(a.n_obs)
lig={"VWF/hemocytin":["113824171"],"THBS4":["113823118"],"POSTN":["113817567","113805211"]}
grp=a.obs.group.astype(str).replace({"blank":"Control","treated":"Hemospasia"})
cts=["MyoEpi","M","MH2","TS"]
fig,axes=plt.subplots(1,3,figsize=(7.2,2.5),sharey=False)
for ax,(name,ids) in zip(axes,lig.items()):
    x=expr(ids); data=[]; pos=[]; labels=[]
    for i,ct in enumerate(cts):
        for j,g in enumerate(["Control","Hemospasia"]):
            m=(a.obs.cell_type==ct).values & (grp==g).values
            data.append(x[m]); pos.append(i+(j-0.5)*0.28); labels.append("")
    bp=ax.boxplot(data,positions=pos,widths=0.24,patch_artist=True,showfliers=False,medianprops=dict(color="#222",lw=0.5))
    for k,p_ in enumerate(bp["boxes"]): p_.set_facecolor("#4C78A8" if k%2==0 else "#E45756"); p_.set_alpha(0.75); p_.set_edgecolor("#444"); p_.set_linewidth(0.3)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts,rotation=30); ax.set_title(name); ax.set_ylabel("log1p CP10K")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
axes[0].legend(handles=[plt.Rectangle((0,0),1,1,color="#4C78A8"),plt.Rectangle((0,0),1,1,color="#E45756")],labels=["Control","Hemospasia"],frameon=False,fontsize=5,loc="upper right")
fig.tight_layout(); save(fig,"Fig4b_ligand_group")

# Fig4 hepatopancreas composition: cluster sizes + PH label-transfer z
g=ad.read_h5ad(os.path.join(OUT,"geo_GSE306674_processed.h5ad"),backed="r")
sizes=g.obs["cl"].astype(str).value_counts().sort_index()
ident={"0":"epithelial/immune","1":"unknown","2":"epithelial/immune","3":"metabolic/absorptive","4":"F cell","5":"metabolic","6":"proliferative","7":"unknown","8":"neural","9":"neural","10":"neural","11":"neural","12":"F cell","13":"F cell","14":"unknown","15":"neuroendocrine","16":"neural"}
colors={"F cell":"#8C510A","metabolic/absorptive":"#41AB5D","metabolic":"#78C679","epithelial/immune":"#B2182B","proliferative":"#E45756","neural":"#2166AC","neuroendocrine":"#7B3294","unknown":"#999999"}
z=pd.read_csv(os.path.join(OUT,"geo_label_transfer_z.csv"),index_col=0)
z.columns=[str(c) for c in z.columns]
phrow=z["2"]
fig,axes=plt.subplots(1,2,figsize=(7.4,3.2),gridspec_kw={"width_ratios":[1,1]})
ax=axes[0]
for i,c in enumerate(sorted(sizes.index,key=int)):
    ax.barh(i,sizes[c],color=colors.get(ident.get(c,"unknown"),"#999"),height=0.62)
    ax.text(sizes[c]+80,i,f"R{c} {ident.get(c,"")}",va="center",fontsize=4.4)
ax.set_yticks([]); ax.invert_yaxis(); ax.set_xlabel("Cells"); ax.set_xlim(0,sizes.max()*1.55)
ax.set_title("GEO hepatopancreas clusters"); 
for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
ax=axes[1]
clusters_sorted=sorted(sizes.index,key=int)
vals=phrow.reindex([f"R{c}" for c in clusters_sorted])
ax.barh(range(len(vals)),vals.values,color=[colors.get(ident.get(c.replace("R",""),"unknown"),"#999") for c in vals.index],height=0.62)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(list(vals.index),fontsize=5); ax.invert_yaxis()
ax.axvline(0,color="#333",lw=0.5); ax.set_xlabel("PH label-transfer z-score")
ax.set_title("PH maps to proliferative R6, not a hemocyte cluster")
for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
fig.tight_layout(); save(fig,"Fig4d_hepatopancreas_composition")
print("done")
