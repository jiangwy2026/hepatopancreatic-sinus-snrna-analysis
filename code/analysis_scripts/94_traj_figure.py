# -*- coding: utf-8 -*-
"""
94_traj_figure.py - pseudotime trajectory figure in the style of the reference panel (D/E):
(D) cells coloured by pseudotime with PH as the root, (E) cells coloured by functional group with
dashed lineage ellipses.

Pseudotime here = shortest-path (geodesic) distance from the PH root in the kNN graph, because the
diffusion map of this dataset has no spectral gap (see logs\traj_rebuild.log), so DPT is uninformative.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
import anndata as ad

OUT = r"D:\桌面\single cell\analysis_output"
FIG = os.path.join(OUT, "figures_pub")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)

ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
LINEAGE = [("PH", "Prohemocyte", "#B2182B"),
           ("GH", "Granulocyte lineage", "#4393C3"),
           ("MH1", "Monocytic hemocyte lineage", "#8C510A")]

def ellipse_from(pts, nsd=2.0):
    mu = pts.mean(axis=0)
    cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov)
    order = ev.argsort()[::-1]
    ev, evec = ev[order], evec[:, order]
    ang = np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), ang

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_fig5c.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"])
    n = A.n_obs
    print("cells:", n, pd.Series(ct).value_counts().to_dict(), flush=True)

    # --- geodesic pseudotime from the PH root
    D = A.obsp["distances"].tocsr() if "distances" in A.obsp else None
    if D is None or D.nnz == 0:
        raise SystemExit("kNN distances missing")
    ph = ct == "PH"
    ph_c = U[ph].mean(axis=0)
    root = int(np.where(ph)[0][np.argmin(((U[ph] - ph_c) ** 2).sum(1))])
    dist = dijkstra(D, directed=False, indices=root)
    fin = np.isfinite(dist)
    print("reachable: %d/%d | raw distance percentiles 50/90/99/max: %s" % (
        int(fin.sum()), n, np.round(np.percentile(dist[fin], [50, 90, 99, 100]), 3).tolist()), flush=True)
    hi = np.percentile(dist[fin], 99.5)
    pt = np.clip(dist / hi, 0, 1)
    pt[~fin] = 1.0
    A.obs["geodesic_pseudotime"] = pt
    tab = pd.DataFrame({"ct": ct, "pt": pt}).groupby("ct")["pt"].agg(["count", "mean", "median"])
    print("--- geodesic pseudotime by group ---"); print(tab.round(3).to_string(), flush=True)

    # --- CytoTRACE-style differentiation potential (independent check)
    ng = np.asarray(A.obs["n_genes_by_counts"], dtype=float)
    if "connectivities" in A.obsp:
        W = A.obsp["connectivities"].tocsr()
        rs = np.asarray(W.sum(1)).ravel(); rs[rs == 0] = 1
        gcs = np.asarray((sp.diags(1 / rs) @ W @ sp.diags(np.ones(n)) * 0).todense()) if False else None
        gcs = np.asarray(W.dot(ng) / rs)
    else:
        gcs = ng
    A.obs["gcs"] = gcs
    print("--- CytoTRACE-style score (smoothed n_genes) by group ---", flush=True)
    print(pd.DataFrame({"ct": ct, "gcs": gcs}).groupby("ct")["gcs"].mean().round(1).to_string(), flush=True)
    A.write_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))

    # --- figure (reference-style D / E)
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.45))
    ax = axes[0]
    o = np.argsort(-pt)
    sc0 = ax.scatter(U[o, 0], U[o, 1], c=pt[o], s=1.1, cmap="Blues_r", linewidths=0,
                     rasterized=True, vmin=0, vmax=1)
    ax.scatter(U[root, 0], U[root, 1], s=26, marker="o", facecolor="none", edgecolor="#B2182B", linewidths=0.9)
    ax.annotate("root (PH)", xy=(U[root, 0], U[root, 1]), xytext=(U[root, 0] + 0.6, U[root, 1] + 2.4),
                fontsize=6, color="#B2182B",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#B2182B"))
    cb = fig.colorbar(sc0, ax=ax, orientation="horizontal", location="top", pad=0.02, fraction=0.05, aspect=32)
    cb.set_label("Pseudotime", fontsize=6.5); cb.ax.tick_params(labelsize=5.5, length=1.4)
    cb.outline.set_visible(False)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    ax = axes[1]
    for g in ORD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.1, c=COL[g], linewidths=0, rasterized=True, label=g)
    for key, lab, col in LINEAGE:
        m = ct == key if key != "MH1" else (ct == "MH1") | (ct == "MH2")
        mu, (w, h), ang = ellipse_from(U[m], 1.9)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (4, 2.2)), lw=0.7, ec=col, zorder=5))
    for key, lab, col in LINEAGE:
        m = ct == key if key != "MH1" else (ct == "MH1") | (ct == "MH2")
        mu = U[m].mean(axis=0)
        ax.text(mu[0], mu[1], lab, fontsize=5.8, color=col, ha="center", va="center", zorder=6,
                path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
    leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=4, frameon=False,
                    markerscale=5, handletextpad=0.3, columnspacing=0.9, fontsize=6)
    for t, g in zip(leg.get_texts(), ORD):
        t.set_text("%s (%d)" % (g, int((ct == g).sum())))
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    for a, letter in zip(axes, ["D.", "E."]):
        a.text(-0.13, 1.08, letter, transform=a.transAxes, fontsize=9, fontweight="bold", va="top")
    for a in axes:
        a.tick_params(labelsize=5.5, length=1.8, width=0.55)
        a.set_xticks([]); a.set_yticks([])
    fig.subplots_adjust(wspace=0.12)
    p = os.path.join(FIG, "Fig4_traj_pseudotime_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight")
    print("saved", p, flush=True)
    json.dump({"root_cell": root, "pseudotime_group_means": tab["mean"].round(4).to_dict(),
               "gcs_group_means": pd.DataFrame({"ct": ct, "gcs": gcs}).groupby("ct")["gcs"].mean().round(2).to_dict()},
              open(os.path.join(OUT, "TRAJ_PSEUDOTIME_REPORT.json"), "w"), indent=1)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_figure.log", "w", encoding="utf8", buffering=1)
    main()
