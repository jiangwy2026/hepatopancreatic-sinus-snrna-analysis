# -*- coding: utf-8 -*-
"""174_pseudotime_smooth_try.py - display-scale DPT/geodesic pseudotime and kNN smoothing."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
from scipy.stats import rankdata

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
FIG = os.path.join(OUT, "figures_pub", "elife_try")
os.makedirs(FIG, exist_ok=True)
MAP = {"MH1": "TS", "PH": "PH", "GH": "GH", "MH2": "MH2"}
ORDER = ["PH", "GH", "TS", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}

def normalize_rank(x):
    x = np.asarray(x, dtype=float)
    good = np.isfinite(x)
    out = np.full(len(x), np.nan)
    if good.sum():
        r = rankdata(x[good], method="average")
        out[good] = (r - 1) / max(1, good.sum() - 1)
    return out

def smooth(x, A, iters=5, alpha=0.35):
    W = A.obsp["connectivities"].tocsr()
    rs = np.asarray(W.sum(1)).ravel(); rs[rs == 0] = 1
    T = sp.diags(1 / rs) @ W
    v = np.asarray(x, dtype=float).copy()
    bad = ~np.isfinite(v)
    v[bad] = np.nanmax(v[~bad]) if (~bad).any() else 1
    for _ in range(iters):
        v = alpha * v + (1 - alpha) * np.asarray(T.dot(v)).ravel()
    return v

def plot(A, vals, title, stem, root):
    U = np.asarray(A.obsm["X_umap"]); ct = A.obs["ct_new"].astype(str).values
    o = np.argsort(-vals)
    fig, ax = plt.subplots(figsize=(3.9, 3.8))
    s = ax.scatter(U[o, 0], U[o, 1], c=vals[o], s=1.3, cmap="viridis", linewidths=0, rasterized=True, vmin=0, vmax=1)
    ax.scatter(U[root, 0], U[root, 1], s=60, marker="*", c="#FF3B30", edgecolors="white", linewidths=0.6, zorder=9)
    ax.annotate("root (PH)", xy=(U[root, 0], U[root, 1]), xytext=(0, 13), textcoords="offset points",
                fontsize=6.0, color="#B2182B", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#B2182B"))
    for g in ORDER:
        c = U[ct == g].mean(0)
        ax.text(c[0], c[1], g, fontsize=6.0, weight="bold", ha="center", va="center", color="white",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="#333333")], zorder=8)
    cb = fig.colorbar(s, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("Pseudotime (display rank)", fontsize=6.2); cb.ax.tick_params(labelsize=5.2); cb.outline.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(title, fontsize=7.2)
    for sp_ in ax.spines.values(): sp_.set_color("#777777"); sp_.set_linewidth(0.5)
    p = os.path.join(FIG, stem + ".png"); fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig); print("saved", p)

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    A.obs["ct_new"] = pd.Categorical(A.obs["ct"].astype(str).map(MAP), categories=ORDER)
    root = int(A.uns.get("iroot", 0))
    ct = A.obs["ct_new"].astype(str).values
    if ct[root] != "PH":
        ph = np.where(ct == "PH")[0]; U = np.asarray(A.obsm["X_umap"])
        root = int(ph[np.argmin(((U[ph] - U[ph].mean(0)) ** 2).sum(1))])
    # DPT rank
    dpt = np.asarray(A.obs["dpt_pseudotime"], dtype=float)
    dpt[dpt == np.inf] = np.nan
    dpt_rank = normalize_rank(dpt)
    plot(A, dpt_rank, "DPT pseudotime (rank display)", "pseudotime_smooth_dpt_rank", root)
    # smoothed DPT
    dpt_s = normalize_rank(smooth(dpt_rank, A, iters=6, alpha=0.30))
    plot(A, dpt_s, "DPT pseudotime, kNN-smoothed", "pseudotime_smooth_dpt", root)
    # geodesic rank and smoothed
    geo = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    geo_rank = normalize_rank(geo)
    plot(A, geo_rank, "Geodesic pseudotime (rank display)", "pseudotime_smooth_geodesic_rank", root)
    geo_s = normalize_rank(smooth(geo_rank, A, iters=6, alpha=0.30))
    plot(A, geo_s, "Geodesic pseudotime, kNN-smoothed", "pseudotime_smooth_geodesic", root)
    # log group medians
    rows=[]
    for name, v in [("dpt",dpt_rank),("dpt_s",dpt_s),("geo",geo_rank),("geo_s",geo_s)]:
        for g in ORDER:
            rows.append({"score":name,"ct":g,"median":float(np.nanmedian(v[ct==g])),"mean":float(np.nanmean(v[ct==g]))})
    print(pd.DataFrame(rows).pivot(index="ct",columns="score",values="median").round(3).to_string())

if __name__ == "__main__":
    log=open(os.path.join(r"D:\桌面\single cell\logs","pseudotime_smooth_try.log"),"w",encoding="utf8",buffering=1)
    sys.stdout=log; sys.stderr=log; main(); log.close()
