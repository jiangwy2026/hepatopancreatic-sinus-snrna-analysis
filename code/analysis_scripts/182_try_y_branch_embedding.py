# -*- coding: utf-8 -*-
"""182_try_y_branch_embedding.py
Exploratory: try to recover a reference-style Y-shaped component layout for
PH -> TS -> GH / MH2. Outputs previews only.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "elife_try")
os.makedirs(PREV, exist_ok=True)
MAP = {"MH1": "TS", "PH": "PH", "GH": "GH", "MH2": "MH2"}
ORDER = ["PH", "GH", "TS", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}

def save(fig, name):
    p = os.path.join(PREV, name)
    fig.savefig(p, dpi=350, bbox_inches="tight", facecolor="white")
    print("saved", p); plt.close(fig)

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    A.obs["ct_new"] = pd.Categorical(A.obs["ct"].astype(str).map(MAP), categories=ORDER)
    ct = A.obs["ct_new"].astype(str).values
    # Balance TS (the largest group) against the smallest informative groups.
    rng = np.random.default_rng(7)
    ts_idx = np.where(ct == "TS")[0]
    take_ts = rng.choice(ts_idx, size=min(len(ts_idx), 900), replace=False)
    keep = np.where(ct != "TS")[0].tolist() + take_ts.tolist()
    keep = np.array(sorted(keep))
    B = A[keep].copy()
    B.obs["ct_new"] = pd.Categorical(B.obs["ct"].astype(str).map(MAP), categories=ORDER)
    print("balanced:", B.n_obs, B.obs["ct_new"].value_counts().to_dict())

    # Balanced diffusion map using the already computed Harmony PCA.
    sc.pp.neighbors(B, use_rep="X_pca_harmony", n_neighbors=20, n_pcs=30)
    sc.tl.diffmap(B, n_comps=15)
    D = np.asarray(B.obsm["X_diffmap"])[:, 1:3]
    bct = B.obs["ct_new"].astype(str).values
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    for g in ORDER:
        m = bct == g
        ax.scatter(D[m, 0], D[m, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True, label=g)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("Component 1", fontsize=6.2); ax.set_ylabel("Component 2", fontsize=6.2)
    ax.legend(frameon=False, fontsize=5.2, markerscale=4, loc="upper right")
    ax.set_title("Balanced diffusion components", fontsize=7)
    for sp in ax.spines.values(): sp.set_color("#777777"); sp.set_linewidth(0.45)
    save(fig, "y_try_balanced_diffmap.png")

    # Try UMAP initialized from a manual Y-shaped PAGA-style layout.
    # Centers: PH lower-left, TS center, GH right, MH2 upper.
    centers = {"PH": np.array([-3.0, -3.5]), "TS": np.array([0.0, 0.0]),
               "GH": np.array([4.2, 0.0]), "MH2": np.array([0.0, 3.5])}
    init = np.zeros((B.n_obs, 2))
    for g in ORDER:
        m = bct == g
        init[m] = centers[g] + rng.normal(0, 0.22, size=(m.sum(), 2))
    sc.tl.umap(B, init_pos=init, random_state=0, min_dist=0.45)
    U = np.asarray(B.obsm["X_umap"])
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    for g in ORDER:
        m = bct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True, label=g)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1", fontsize=6.2); ax.set_ylabel("UMAP 2", fontsize=6.2)
    ax.legend(frameon=False, fontsize=5.2, markerscale=4, loc="upper right")
    ax.set_title("Y-init UMAP (balanced)", fontsize=7)
    for sp in ax.spines.values(): sp.set_color("#777777"); sp.set_linewidth(0.45)
    save(fig, "y_try_init_umap.png")

if __name__ == "__main__":
    log = open(r"D:\桌面\single cell\logs\y_branch_try.log", "w", encoding="utf8", buffering=1)
    sys.stdout = log; sys.stderr = log; main(); log.close()
