# -*- coding: utf-8 -*-
"""248_step2_trajectory_projection.py

Step 2: internal Harmony integration of hepatic-sinusoid niche + hemocytes +
hepatopancreas, then project all cells onto the hemocyte differentiation
trajectory (PH -> GH) and quantify progenitor proximity.

Outputs in analysis_output/step2/.
"""
import os, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import harmonypy as hm
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\桌面\single cell\analysis_output"
STEP = os.path.join(OUT, "step2")
os.makedirs(STEP, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT, ".mplconfig"))

GROUP = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
NICHE = ["TS","MH2","MyoEpi","M"]
HEMO = ["PH","GH"]
HEPATO = ["F","R","B"]


def main():
    a = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    a.X = a.layers["counts"].copy()
    a.obs["cell_type"] = a.obs["vendor_cluster"].astype(str).map(GROUP)
    a.obs["compartment"] = a.obs["cell_type"].map(
        lambda x: "niche" if x in NICHE else ("hemocyte" if x in HEMO else ("hepatopancreas" if x in HEPATO else "other")))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a.layers["lognorm"] = a.X.copy()

    # HVG + PCA + Harmony (batch = sample)
    b = a.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, batch_key="sample")
    b = b[:, b.var["highly_variable"].values].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=50, svd_solver="arpack", random_state=0)
    ho = hm.run_harmony(b.obsm["X_pca"], b.obs, ["sample"], random_state=0, max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr).T
    a.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=20, random_state=0)
    sc.tl.umap(a, random_state=0)

    # hemocyte trajectory: PH -> GH
    h = a[a.obs["cell_type"].isin(HEMO)].copy()
    sc.pp.neighbors(h, use_rep="X_pca_harmony", n_neighbors=20, random_state=0)
    sc.tl.diffmap(h, n_comps=15)
    # root: PH cell with highest ccnd2
    sym = h.var["sym"].astype(str).str.lower()
    cc_vars = h.var_names[sym.isin(["ccnd2"]).values].tolist()
    if cc_vars:
        v = np.asarray(h[:, cc_vars].X.todense()).mean(axis=1)
    else:
        v = np.zeros(h.n_obs)
    ph = (h.obs["cell_type"] == "PH").values
    root_idx = int(np.where(ph)[0][np.argmax(v[ph])])
    root_ccnd2 = float(v[ph].max())
    h.uns["iroot"] = root_idx
    sc.tl.dpt(h, n_dcs=10)
    a.obs["dpt_hemocyte"] = np.nan
    a.obs.loc[h.obs_names, "dpt_hemocyte"] = h.obs["dpt_pseudotime"].values

    # project non-hemocyte cells to hemocyte manifold
    idx_h = np.where(a.obs["cell_type"].isin(HEMO).values)[0]
    Z_h = a.obsm["X_pca_harmony"][idx_h]
    nn = NearestNeighbors(n_neighbors=30).fit(Z_h)
    d, ind = nn.kneighbors(a.obsm["X_pca_harmony"])
    pt = h.obs["dpt_pseudotime"].values[ind]
    w = 1.0 / (d + 1e-6)
    proj = (pt * w).sum(1) / w.sum(1)
    a.obs["projected_pseudotime"] = proj
    a.obs["proj_pseudotime"] = np.where(a.obs["cell_type"].isin(HEMO).values,
                                        a.obs["dpt_hemocyte"].values, proj)
    cell_arr = np.asarray(h.obs["cell_type"].astype(str))
    ph_neighbor = (cell_arr[ind] == "PH").mean(1)
    gh_neighbor = (cell_arr[ind] == "GH").mean(1)
    a.obs["PH_neighbor_frac"] = ph_neighbor
    a.obs["GH_neighbor_frac"] = gh_neighbor

    # receptor / niche reception score (SonicParanoid-confirmed receptors)
    rec_gids = ["113830322","113805715","113800203","113821319","113827435","113802681"]
    g2v = {str(g).replace("ncbi_",""): g for g in a.var_names}
    rec_vars = [g2v[g] for g in rec_gids if g in g2v]
    Xr = np.asarray(a[:, rec_vars].layers["lognorm"].todense(), dtype=float)
    Zr = (Xr - Xr.mean(0)) / (Xr.std(0) + 1e-9)
    a.obs["receptor_score"] = Zr.mean(1)

    a.write_h5ad(os.path.join(STEP, "step2_integrated.h5ad"), compression="gzip")
    cols = ["sample","cell_type","compartment","dpt_hemocyte","projected_pseudotime",
            "proj_pseudotime","PH_neighbor_frac","GH_neighbor_frac","receptor_score"]
    a.obs[cols].to_csv(os.path.join(STEP, "step2_cell_scores.csv"), encoding="utf-8-sig")

    # summary by cell type
    rows = []
    for ct, g in a.obs.groupby("cell_type"):
        rows.append({"cell_type":ct,"n":len(g),
                     "median_projected_pseudotime":float(g["projected_pseudotime"].median()),
                     "median_PH_neighbor_frac":float(g["PH_neighbor_frac"].median()),
                     "median_GH_neighbor_frac":float(g["GH_neighbor_frac"].median()),
                     "median_receptor_score":float(g["receptor_score"].median())})
    pd.DataFrame(rows).sort_values("median_projected_pseudotime").to_csv(
        os.path.join(STEP,"step2_summary_by_celltype.csv"), index=False, encoding="utf-8-sig")

    # figures
    fig, axes = plt.subplots(1, 3, figsize=(15,4.2))
    for ax, key, title in zip(axes, ["cell_type","compartment","proj_pseudotime"],
                              ["Cell type","Compartment","Projected pseudotime"]):
        if key == "proj_pseudotime":
            v = a.obs[key].values
            sca = ax.scatter(a.obsm["X_umap"][:,0], a.obsm["X_umap"][:,1], c=v, cmap="viridis", s=3, rasterized=True)
            plt.colorbar(sca, ax=ax, fraction=0.04)
        else:
            cat = a.obs[key].astype("category")
            ax.scatter(a.obsm["X_umap"][:,0], a.obsm["X_umap"][:,1], c=cat.cat.codes, cmap="tab20", s=3, rasterized=True)
            for code, name in enumerate(cat.cat.categories):
                m = cat.cat.codes == code
                if m.sum() > 20:
                    ax.text(np.median(a.obsm["X_umap"][m.values,0]), np.median(a.obsm["X_umap"][m.values,1]),
                            name, fontsize=7, ha="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step2_integrated_umap.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step2_integrated_umap.pdf"),bbox_inches="tight"); plt.close(fig)

    order = ["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
    fig, axes = plt.subplots(1, 2, figsize=(11,4.2))
    data = [a.obs.loc[a.obs.cell_type==c,"projected_pseudotime"].dropna().values for c in order]
    axes[0].boxplot(data, tick_labels=order, showfliers=False)
    axes[0].set_title("Projected pseudotime (hemocyte trajectory)", fontsize=9)
    axes[0].tick_params(axis="x", rotation=45)
    data2 = [a.obs.loc[a.obs.cell_type==c,"receptor_score"].dropna().values for c in order]
    axes[1].boxplot(data2, tick_labels=order, showfliers=False)
    axes[1].set_title("Niche-receptor score (ITGB1/ITGAV/ITGA4/NOTCH1)", fontsize=9)
    axes[1].tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step2_pseudotime_receptor_boxplot.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step2_pseudotime_receptor_boxplot.pdf"),bbox_inches="tight"); plt.close(fig)

    summary = {"n_cells": int(a.n_obs),
               "groups": {k:int(v) for k,v in a.obs.cell_type.value_counts().items()},
               "root_cell_type": "PH", "root_ccnd2": root_ccnd2}
    with open(os.path.join(STEP,"step2_summary.json"),"w",encoding="utf-8") as f:
        json.dump(summary,f,ensure_ascii=False,indent=2)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
