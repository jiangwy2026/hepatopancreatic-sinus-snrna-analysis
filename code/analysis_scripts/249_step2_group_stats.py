import numpy as np, pandas as pd, anndata as ad, json, sys, os
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")
OUT = r"D:\桌面\single cell\analysis_output"; STEP = os.path.join(OUT,"step2")
a = ad.read_h5ad(os.path.join(STEP,"step2_integrated.h5ad"))
rec_gids = ["113830322","113805715","113800203","113821319","113827435","113802681"]
g2v = {str(g).replace("ncbi_",""): g for g in a.var_names}
rec_vars = [g2v[g] for g in rec_gids if g in g2v]
X = np.asarray(a[:, rec_vars].layers["lognorm"].todense(), dtype=float)
pos = (X > 0).any(axis=1)
order = ["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
rows=[]
for ct in order:
    m = (a.obs.cell_type == ct).values
    rows.append({"cell_type":ct,"n":int(m.sum()),
        "mean_pseudotime":float(a.obs.projected_pseudotime[m].mean()),
        "median_pseudotime":float(a.obs.projected_pseudotime[m].median()),
        "mean_PH_neighbor_frac":float(a.obs.PH_neighbor_frac[m].mean()),
        "frac_receptor_positive":float(pos[m].mean()),
        "mean_ITGB1":float(X[m,0].mean()),"mean_ITGAV":float(X[m,1].mean()),
        "mean_ITGA4":float(X[m,2].mean()),"mean_NOTCH1":float(X[m,3].mean())})
df=pd.DataFrame(rows); df.to_csv(os.path.join(STEP,"step2_group_stats.csv"),index=False,encoding="utf-8-sig")
# tests
kw_pt = stats.kruskal(*[a.obs.projected_pseudotime[(a.obs.cell_type==c).values].values for c in order])
kw_ph = stats.kruskal(*[a.obs.PH_neighbor_frac[(a.obs.cell_type==c).values].values for c in order])
pair=[]
for c in order:
    if c=="PH": continue
    u,p = stats.mannwhitneyu(a.obs.PH_neighbor_frac[(a.obs.cell_type=="PH").values],
                             a.obs.PH_neighbor_frac[(a.obs.cell_type==c).values], alternative="two-sided")
    pair.append({"group":c,"PH_vs_group_p":float(p)})
tests={"kruskal_pseudotime_p":float(kw_pt.pvalue),"kruskal_PH_neighbor_p":float(kw_ph.pvalue),"pairwise":pair}
json.dump(tests,open(os.path.join(STEP,"step2_tests.json"),"w",encoding="utf-8"),ensure_ascii=False,indent=2)
pd.set_option("display.width",250); print(df.to_string(index=False)); print(tests)
# plots
fig,axes=plt.subplots(1,3,figsize=(15,4.2))
axes[0].boxplot([a.obs.PH_neighbor_frac[(a.obs.cell_type==c).values] for c in order],tick_labels=order,showfliers=False)
axes[0].set_title("Fraction of PH among 30 nearest hemocytes",fontsize=9); axes[0].tick_params(axis="x",rotation=45)
axes[1].bar(order,[df.set_index("cell_type").loc[c,"frac_receptor_positive"] for c in order],color="#C44E52")
axes[1].set_title("Fraction receptor-positive (ITGB1/ITGAV/ITGA4/NOTCH1)",fontsize=9); axes[1].tick_params(axis="x",rotation=45)
M=df.set_index("cell_type").loc[order,["mean_ITGB1","mean_ITGAV","mean_ITGA4","mean_NOTCH1"]].values
im=axes[2].imshow(M.T,aspect="auto",cmap="Reds")
axes[2].set_xticks(range(len(order))); axes[2].set_xticklabels(order,rotation=45)
axes[2].set_yticks(range(4)); axes[2].set_yticklabels(["ITGB1","ITGAV","ITGA4","NOTCH1"])
axes[2].set_title("Mean receptor expression (log1p CP10K)",fontsize=9); plt.colorbar(im,ax=axes[2],fraction=0.045)
fig.tight_layout(); fig.savefig(os.path.join(STEP,"step2_group_comparison.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(STEP,"step2_group_comparison.pdf"),bbox_inches="tight"); plt.close(fig)
