# -*- coding: utf-8 -*-
"""173_pseudotime_try_elife.py
Try several reference-style pseudotime layouts on the four hemocyte states.
Outputs are exploratory previews only.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

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


def save(fig, stem, dpi=300):
    p = os.path.join(FIG, stem + ".png")
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    print("saved", p)
    plt.close(fig)


def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).map(MAP).values
    A.obs["ct_new"] = pd.Categorical(ct, categories=ORDER)
    U = np.asarray(A.obsm["X_umap"])
    DC = np.asarray(A.obsm["X_diffmap"])[:, 1:4]
    root = int(A.uns.get("iroot", 0))
    # root should be in PH; if not, use the PH cell closest to the PH centroid
    if ct[root] != "PH":
        ph = np.where(ct == "PH")[0]
        root = int(ph[np.argmin(((U[ph] - U[ph].mean(0)) ** 2).sum(1))])
    A.uns["iroot"] = root
    print("root", root, A.obs_names[root], "ct", ct[root])
    print("dpt before", np.isfinite(np.asarray(A.obs["dpt_pseudotime"], dtype=float)).sum())
    # recompute DPT on the existing diffusion map
    sc.tl.dpt(A, n_dcs=10)
    pt = np.asarray(A.obs["dpt_pseudotime"], dtype=float)
    good = np.isfinite(pt)
    print("dpt finite", good.sum(), "range", np.nanmin(pt), np.nanmax(pt))
    print(pd.DataFrame({"ct": ct, "pt": pt}).groupby("ct")["pt"].agg(["count", "median", "mean"]).round(3).to_string())
    A.obs["dpt_recomputed"] = pt

    # shared axes limits
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.04 * (x1 - x0); py = 0.04 * (y1 - y0)
    d0, d1 = np.percentile(DC[:, 0], [0.3, 99.7]); e0, e1 = np.percentile(DC[:, 1], [0.3, 99.7])

    def style(ax, xlabel, ylabel):
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(xlabel, fontsize=6.5); ax.set_ylabel(ylabel, fontsize=6.5)
        for sp in ax.spines.values():
            sp.set_color("#777777"); sp.set_linewidth(0.5)

    # 1) diffusion component plot colored by DPT
    fig, ax = plt.subplots(figsize=(3.9, 3.8))
    o = np.argsort(-pt)
    s = ax.scatter(DC[o, 0], DC[o, 1], c=pt[o], s=1.4, cmap="viridis", linewidths=0,
                   rasterized=True, vmin=np.nanpercentile(pt, 1), vmax=np.nanpercentile(pt, 99))
    ax.scatter(DC[root, 0], DC[root, 1], s=55, marker="*", c="#E31A1C", edgecolors="white", linewidths=0.6, zorder=7)
    ax.annotate("root (PH)", xy=(DC[root, 0], DC[root, 1]), xytext=(0, 12), textcoords="offset points",
                fontsize=6, color="#B2182B", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#B2182B"))
    cb = fig.colorbar(s, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("DPT pseudotime", fontsize=6.2); cb.ax.tick_params(labelsize=5.2); cb.outline.set_visible(False)
    style(ax, "Diffusion component 1", "Diffusion component 2")
    ax.set_title("Diffusion map + PH-rooted DPT", fontsize=7.5)
    save(fig, "pseudotime_try_diffmap")

    # 2) UMAP colored by DPT
    fig, ax = plt.subplots(figsize=(3.9, 3.8))
    s = ax.scatter(U[o, 0], U[o, 1], c=pt[o], s=1.4, cmap="viridis", linewidths=0,
                   rasterized=True, vmin=np.nanpercentile(pt, 1), vmax=np.nanpercentile(pt, 99))
    ax.scatter(U[root, 0], U[root, 1], s=55, marker="*", c="#E31A1C", edgecolors="white", linewidths=0.6, zorder=7)
    ax.annotate("root (PH)", xy=(U[root, 0], U[root, 1]), xytext=(0, 12), textcoords="offset points",
                fontsize=6, color="#B2182B", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#B2182B"))
    cb = fig.colorbar(s, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("DPT pseudotime", fontsize=6.2); cb.ax.tick_params(labelsize=5.2); cb.outline.set_visible(False)
    style(ax, "UMAP 1", "UMAP 2")
    ax.set_title("UMAP + PH-rooted DPT", fontsize=7.5)
    save(fig, "pseudotime_try_umap")

    # 3) UMAP colored by group, with group centroids and PAGA edges
    fig, ax = plt.subplots(figsize=(3.9, 3.8))
    for g in ORDER:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.4, c=COL[g], linewidths=0, rasterized=True, label="%s (%d)" % (g, int(m.sum())))
    cent = {g: U[ct == g].mean(0) for g in ORDER}
    # PAGA connectivity as a descriptive summary; no arrow direction is inferred here
    sc.tl.paga(A, groups="ct_new")
    C = pd.DataFrame(A.uns["paga"]["connectivities"].toarray(), index=A.obs["ct_new"].cat.categories,
                     columns=A.obs["ct_new"].cat.categories)
    print("PAGA connectivity\n", C.round(3).to_string())
    for i in range(len(ORDER)):
        for j in range(i + 1, len(ORDER)):
            w = float(C.iloc[i, j])
            if w > 0.05:
                p1, p2 = cent[ORDER[i]], cent[ORDER[j]]
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#666666", lw=0.4 + 4.0 * w,
                        alpha=0.45, zorder=4, solid_capstyle="round")
    for g in ORDER:
        c = cent[g]
        ax.text(c[0], c[1], g, fontsize=6.2, weight="bold", ha="center", va="center", color="black",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")], zorder=8)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.02), ncol=2, frameon=False,
              markerscale=4.0, fontsize=5.2, handletextpad=0.25, columnspacing=0.7)
    style(ax, "UMAP 1", "UMAP 2")
    ax.set_title("Groups + PAGA connectivity", fontsize=7.5)
    save(fig, "pseudotime_try_paga")

    # 4) UMAP + geodesic pseudotime + PH-rooted path through group centroids
    raw = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    # normalized rank for display
    rankpt = np.argsort(np.argsort(raw)) / (len(raw) - 1.0)
    fig, ax = plt.subplots(figsize=(3.9, 3.8))
    s = ax.scatter(U[o, 0], U[o, 1], c=rankpt[o], s=1.4, cmap="magma", linewidths=0,
                   rasterized=True, vmin=0, vmax=1)
    ax.scatter(U[root, 0], U[root, 1], s=55, marker="*", c="#00FFFF", edgecolors="black", linewidths=0.55, zorder=8)
    for g in ORDER:
        c = U[ct == g].mean(0)
        ax.text(c[0], c[1], g, fontsize=6.0, weight="bold", ha="center", va="center", color="white",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="#333333")], zorder=9)
    cb = fig.colorbar(s, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("Pseudotime (rank)", fontsize=6.2); cb.ax.tick_params(labelsize=5.2); cb.outline.set_visible(False)
    style(ax, "UMAP 1", "UMAP 2")
    ax.set_title("Geodesic pseudotime from PH", fontsize=7.5)
    save(fig, "pseudotime_try_geodesic")

    # 5) combined 2x2 exploratory board
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.55))
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1" if ax is not axes[0] else "Diffusion component 1")
        ax.set_ylabel("UMAP 2" if ax is not axes[0] else "Diffusion component 2")
    # re-draw compactly
    ax = axes[0]
    ax.scatter(DC[o, 0], DC[o, 1], c=pt[o], s=0.8, cmap="viridis", linewidths=0, rasterized=True)
    ax.scatter(DC[root, 0], DC[root, 1], s=28, marker="*", c="#E31A1C", edgecolors="white", linewidths=0.4, zorder=7)
    ax.set_title("Diffusion map + DPT", fontsize=6.6)
    ax = axes[1]
    ax.scatter(U[o, 0], U[o, 1], c=pt[o], s=0.8, cmap="viridis", linewidths=0, rasterized=True)
    ax.scatter(U[root, 0], U[root, 1], s=28, marker="*", c="#E31A1C", edgecolors="white", linewidths=0.4, zorder=7)
    ax.set_title("UMAP + DPT", fontsize=6.6)
    ax = axes[2]
    for g in ORDER:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=0.8, c=COL[g], linewidths=0, rasterized=True)
    for i in range(len(ORDER)):
        for j in range(i + 1, len(ORDER)):
            w = float(C.iloc[i, j])
            if w > 0.05:
                p1, p2 = cent[ORDER[i]], cent[ORDER[j]]
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#666666", lw=0.3 + 3.0 * w, alpha=0.5, zorder=4)
    for g in ORDER:
        c = cent[g]
        ax.text(c[0], c[1], g, fontsize=5.2, weight="bold", ha="center", va="center",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=8)
    ax.set_title("Groups + PAGA", fontsize=6.6)
    fig.tight_layout(w_pad=0.5)
    save(fig, "pseudotime_try_board")


if __name__ == "__main__":
    log = open(os.path.join(r"D:\桌面\single cell\logs", "pseudotime_try.log"), "w", encoding="utf8", buffering=1)
    sys.stdout = log
    sys.stderr = log
    main()
    log.close()
