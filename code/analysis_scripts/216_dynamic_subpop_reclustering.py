# -*- coding: utf-8 -*-
"""216_dynamic_subpop_reclustering.py
Recluster the five composition-changing subpopulations, save subcluster UMAPs,
subcluster composition, and top markers. Preview-only figures.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
SUB_DIR = os.path.join(OUT, "dynamic_subpop_subclusters")
os.makedirs(PARTS, exist_ok=True); os.makedirs(SUB_DIR, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style

apply_pub_style(7)
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
sc.settings.verbosity = 0

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
DYN = [("PH", 2.30), ("MH2", 2.17), ("R", 2.29), ("MyoEpi", 0.44), ("PM-like", 0.27)]
PALETTE = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#9C755F"]


def recluster_one(A, grp, fold):
    sub = A[A.obs["ct"].astype(str) == grp].copy()
    sub.obs["group2"] = sub.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    sc.pp.filter_genes(sub, min_cells=3)
    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.highly_variable_genes(sub, n_top_genes=min(1500, sub.n_vars), flavor="seurat")
    sub = sub[:, sub.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    n_pcs = min(20, max(5, sub.n_obs // 15))
    sc.tl.pca(sub, n_comps=n_pcs, svd_solver="arpack", random_state=0)
    # Integrate Control and Hemospasia before subclustering so the dominant axis
    # is cell state rather than sample group.
    sc.external.pp.harmony_integrate(sub, key="group2", basis="X_pca",
                                     adjusted_basis="X_pca_harmony", random_state=0)
    sc.pp.neighbors(sub, n_neighbors=min(15, max(5, sub.n_obs // 20)),
                    use_rep="X_pca_harmony", random_state=0)
    res = {"PH": 0.50, "MH2": 0.60, "R": 0.35,
           "MyoEpi": 0.30, "PM-like": 0.03}[grp]
    sc.tl.leiden(sub, resolution=res, random_state=0, flavor="igraph", n_iterations=2,
                 directed=False, key_added="subcluster")
    ncl = sub.obs["subcluster"].nunique()
    if ncl < 2:
        sc.tl.leiden(sub, resolution=0.45, random_state=0, flavor="igraph", n_iterations=2,
                     directed=False, key_added="subcluster")
    elif ncl > 8:
        sc.tl.leiden(sub, resolution=0.05, random_state=0, flavor="igraph", n_iterations=2,
                     directed=False, key_added="subcluster")
    sc.tl.umap(sub, random_state=0)
    sub.obs["subcluster"] = pd.Categorical(sub.obs["subcluster"].astype(str),
                                           categories=sorted(sub.obs["subcluster"].astype(str).unique(),
                                                             key=lambda x: int(x)))
    sc.tl.rank_genes_groups(sub, "subcluster", method="wilcoxon", tie_correct=False)
    sym = dict(zip(sub.var_names.astype(str), sub.var["sym"].astype(str)))
    rows = []
    for cl in sub.obs["subcluster"].cat.categories:
        d = sc.get.rank_genes_groups_df(sub, group=str(cl)).copy()
        d["gid"] = d["names"].astype(str)
        d["sym"] = d["gid"].map(sym)
        d = d[np.isfinite(d.pvals_adj)].sort_values(["pvals_adj", "scores"], ascending=[True, False]).head(10)
        for r in d.itertuples():
            rows.append({"cluster": grp, "subcluster": str(cl), "gene": r.sym, "gid": r.gid,
                         "score": float(r.scores), "log2FC": float(r.logfoldchanges),
                         "pvals_adj": float(r.pvals_adj)})
    pd.DataFrame(rows).to_csv(os.path.join(SUB_DIR, f"{grp}_subcluster_markers.csv"),
                              index=False, encoding="utf-8-sig")
    obs = sub.obs[["group2", "subcluster"]].copy()
    obs["umap1"] = sub.obsm["X_umap"][:, 0]
    obs["umap2"] = sub.obsm["X_umap"][:, 1]
    obs.to_csv(os.path.join(SUB_DIR, f"{grp}_subcluster_umap.csv"), encoding="utf-8-sig")
    return sub, rows


def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.obs["ct"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    results = {}
    for grp, fold in DYN:
        print("Reclustering", grp)
        sub, rows = recluster_one(A, grp, fold)
        results[grp] = sub
        print(f"  {grp}: {sub.n_obs} cells, subclusters {sub.obs['subcluster'].nunique()}, "
              f"counts {sub.obs['subcluster'].astype(str).value_counts().to_dict()}")
    comp_rows = []
    for grp, fold in DYN:
        sub = results[grp]
        tab = pd.crosstab(sub.obs["subcluster"].astype(str), sub.obs["group2"].astype(str))
        tab = tab.reindex(columns=["Control", "Hemospasia"], fill_value=0)
        for sc_id, r in tab.iterrows():
            nc, nh = int(r["Control"]), int(r["Hemospasia"])
            pc = 100 * nc / tab["Control"].sum() if tab["Control"].sum() else 0
            ph = 100 * nh / tab["Hemospasia"].sum() if tab["Hemospasia"].sum() else 0
            comp_rows.append({"population": grp, "subcluster": str(sc_id),
                              "n_control": nc, "n_hemospasia": nh,
                              "pct_control": pc, "pct_hemospasia": ph,
                              "fold_hemospasia_vs_control": ph / pc if pc > 0 else np.nan})
    pd.DataFrame(comp_rows).to_csv(os.path.join(OUT, "Fig5_dynamic_subpop_subcluster_composition.csv"),
                                   index=False, encoding="utf-8-sig")

    # UMAP preview
    fig = plt.figure(figsize=(7.2, 5.3), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=0.05, right=0.97, top=0.90, bottom=0.08,
                          wspace=0.22, hspace=0.35)
    pos = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    for (grp, fold), p in zip(DYN, pos):
        ax = fig.add_subplot(gs[p])
        sub = results[grp]
        U = sub.obsm["X_umap"]
        cl = sub.obs["subcluster"].astype(str).values
        cats = list(sub.obs["subcluster"].cat.categories)
        for k, c in enumerate(cats):
            m = cl == c
            ax.scatter(U[m, 0], U[m, 1], s=2.2, c=PALETTE[k % len(PALETTE)],
                       linewidths=0, rasterized=True, alpha=0.85)
            if m.sum() >= 3:
                cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
                ax.text(cx, cy, c, fontsize=4.8, ha="center", va="center",
                        color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.12", fc=PALETTE[k % len(PALETTE)],
                                  ec="none", alpha=0.75))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.35); sp.set_color("#777777")
        ax.set_title(f"{grp}  {'↑' if fold > 1 else '↓'} {fold:.2f}×  ({sub.n_obs} cells)",
                     fontsize=5.8, pad=3, color="#B2182B" if fold > 1 else "#2166AC",
                     fontweight="bold")
    fig.text(0.5, 0.965, "Subclustering of the five dynamic populations", ha="center",
             va="top", fontsize=7.2)
    out1 = os.path.join(PARTS, "Fig5_dynamic_subpop_subclusters_UMAP_preview.png")
    fig.savefig(out1, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out1)
    plt.close(fig)

    # Composition preview: two stacked bars per population
    fig, axes = plt.subplots(1, 5, figsize=(7.2, 2.75), facecolor="white")
    for ax, (grp, fold) in zip(axes, DYN):
        sub = results[grp]
        tab = pd.crosstab(sub.obs["subcluster"].astype(str), sub.obs["group2"].astype(str))
        tab = tab.reindex(columns=["Control", "Hemospasia"], fill_value=0)
        pct = tab.div(tab.sum(axis=0), axis=1).fillna(0) * 100
        cats = list(pct.index)
        left = np.zeros(2)
        for k, c in enumerate(cats):
            v = pct.loc[c, ["Control", "Hemospasia"]].values
            ax.barh([1, 0], v, left=left, color=PALETTE[k % len(PALETTE)],
                    height=0.55, edgecolor="white", linewidth=0.4)
            for j, vv in enumerate(v):
                if vv >= 12:
                    ax.text(left[j] + vv / 2, 1 - j, f"{vv:.0f}", ha="center", va="center",
                            fontsize=4.5, color="white")
            left += v
        ax.set_xlim(0, 100); ax.set_ylim(-0.55, 1.55)
        ax.set_yticks([1, 0]); ax.set_yticklabels(["Ctrl", "HSP"], fontsize=5.0)
        ax.set_xticks([0, 50, 100]); ax.tick_params(axis="x", labelsize=4.8, length=1.2, width=0.4)
        ax.tick_params(axis="y", length=0, pad=1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_title(f"{grp}\n{'↑' if fold > 1 else '↓'} {fold:.2f}×", fontsize=5.5,
                     color="#B2182B" if fold > 1 else "#2166AC", fontweight="bold", pad=2)
    max_n = max(results[g].obs["subcluster"].nunique() for g, _ in DYN)
    handles = [Patch(color=PALETTE[k % len(PALETTE)], label=f"s{k}") for k in range(max_n)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=min(6, max_n),
               frameon=False, fontsize=5.0, handlelength=0.7, handletextpad=0.3,
               columnspacing=1.0)
    fig.suptitle("Subcluster composition by group", y=1.02, fontsize=7.0)
    fig.text(0.5, -0.19, "Harmony-integrated by group before subclustering. Bars show subcluster percentages within Control (Ctrl) and Hemospasia (HSP).",
             ha="center", va="top", fontsize=5.0, color="#555555")
    out2 = os.path.join(PARTS, "Fig5_dynamic_subpop_subclusters_composition_preview.png")
    fig.savefig(out2, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out2)
    plt.close(fig)


if __name__ == "__main__":
    main()
