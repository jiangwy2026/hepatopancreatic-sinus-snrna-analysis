# -*- coding: utf-8 -*-
"""168_blood_figs_rename.py - 血细胞气泡图、细胞周期、谱系图（MH1→TS 重出）"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from scipy.interpolate import UnivariateSpline
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "TS", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}
GROUPS = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}

def dotplot():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2":"PH","1":"GH","0":"TS","3":"MH2"}
    ct = A.obs["vendor_cluster"].astype(str).map(base)
    H = A[ct.isin(ORD).values].copy()
    H.obs["ct"] = pd.Categorical(ct[ct.isin(ORD)].values, categories=ORD)
    H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    gpos = {g: i for i, g in enumerate(H.var_names.astype(str))}
    P = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); P["id"] = P["gid"].astype(str)
    X = H[:, [gpos[g] for g in P.id]].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    pct = np.zeros((len(P), 4)); mean = np.zeros_like(pct)
    for j, g in enumerate(ORD):
        m = (H.obs["ct"] == g).values
        pct[:, j] = 100 * (X[m] > 0).mean(0); mean[:, j] = X[m].mean(0)
    z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)
    blocks = [("Prohemocyte", P[P.grp == "PH"]), ("Granulocyte", P[P.grp == "GH"]),
              ("Monocytic / transition", P[P.grp.isin(["MH1", "MH2"])])]
    order = np.array([i for blk in blocks for i in blk[1].index])
    labels, pcto, zo = P.sym.values[order], pct[order], z[order]
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    for i in range(len(labels)):
        for j in range(len(ORD)):
            if pcto[i, j] <= 0.5: continue
            ax.scatter(i, j, s=3 + (pcto[i, j] / 100) ** 1.25 * 62,
                       c=[plt.cm.RdBu_r((np.clip(zo[i, j], -2, 2) + 2) / 4)],
                       edgecolors="#555555", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.7, len(labels) - 0.3); ax.set_ylim(len(ORD) - 0.45, -0.55)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.2, style="italic")
    ax.set_yticks(range(len(ORD))); ax.set_yticklabels(ORD, fontsize=7, fontweight="bold")
    ax.tick_params(axis="both", length=0, pad=2)
    ax.grid(True, ls=(0, (1, 2)), lw=0.35, color="#CCCCCC", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)
    k = 0
    for name, blk in blocks:
        a, b = k, k + len(blk) - 1; k = b + 1
        ax.plot([a - 0.45, b + 0.45], [-0.72, -0.72], color="black", lw=0.9, clip_on=False, zorder=4)
        ax.text((a + b) / 2, -1.05, name, ha="center", va="bottom", fontsize=6.2, weight="bold")
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + 0.42 * pos.height, 0.010, 0.34 * pos.height])
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-2, 2))
    cb = fig.colorbar(sm, cax=cax); cb.set_label("avg.exp.scale", fontsize=5.6, labelpad=1)
    cb.ax.tick_params(labelsize=5, length=1.2); cb.outline.set_visible(False)
    sx = pos.x1 + 0.055
    for kk, p in enumerate((20, 40, 60, 80)):
        y = 0.95 - kk * 0.16
        fig.text(sx, pos.y0 + y * pos.height, "●", fontsize=4 + (p / 100) ** 1.25 * 62 ** 0.5, va="center", color="#666666")
        fig.text(sx + 0.030, pos.y0 + y * pos.height, "%d" % p, fontsize=5.2, va="center")
    fig.text(sx + 0.010, pos.y0 + 1.02 * pos.height, "pct.exp", fontsize=5.6, va="bottom")
    p = os.path.join(FIG, "Fig_blood4_dotplot_TS_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

def cellcycle():
    comp = pd.read_csv(os.path.join(OUT, "cellcycle_phase_composition.csv"), index_col=0)
    comp = comp.rename(index={"MH1": "TS"}).reindex(ORD)
    PHCOL = {"G1": "#BFBFBF", "S": "#4393C3", "G2M": "#B2182B"}
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    left = np.zeros(len(ORD))
    for ph in ("G1", "S", "G2M"):
        v = comp[ph].values
        ax.barh(np.arange(len(ORD)), v, left=left, color=PHCOL[ph], edgecolor="white", height=0.62, label=ph)
        for k, (l, x) in enumerate(zip(left, v)):
            if x >= 5: ax.text(l + x / 2, k, "%.0f" % x, fontsize=6, ha="center", va="center",
                                color="#333333" if ph == "G1" else "white")
        left += v
    ax.set_xlim(0, 100); ax.set_ylim(len(ORD) - 0.45, -0.55)
    ax.set_yticks(range(len(ORD))); ax.set_yticklabels(ORD, fontsize=8, fontweight="bold")
    ax.set_xlabel("Cells per phase (%)", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=6.5, length=1.8, width=0.5); ax.tick_params(axis="y", length=0, pad=3)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False, fontsize=6.5,
              handlelength=0.9, handletextpad=0.3, columnspacing=0.9)
    p = os.path.join(FIG, "Fig_cellcycle_bar_TS_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def lineage(U, ct, out, xlabel, title, legend=True):
    fig, ax = plt.subplots(figsize=(3.7, 4.0))
    for g in ORD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.3, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())))
    for m, off in [((ct == "GH"), (0.4, 4.3)), ((ct == "MH2") | (ct == "TS"), (-5.4, -5.9)), ((ct == "PH"), (1.4, 4.6))]:
        mu, (w, h), ang = ell(U[m], 1.7)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
    cent = {g: U[ct == g].mean(0) for g in ORD}
    for path, c in [(("PH", "GH"), "#2E6DA4"), (("PH", "MH2", "TS"), "#7F4A12")]:
        p = np.vstack([cent[g] for g in path])
        t = np.arange(len(p), dtype=float)
        fx = UnivariateSpline(t, p[:, 0], k=min(3, len(p) - 1)); fy = UnivariateSpline(t, p[:, 1], k=min(3, len(p) - 1))
        tt = np.linspace(0, len(p) - 1, 220); curve = np.c_[fx(tt), fy(tt)]
        ax.plot(curve[:, 0], curve[:, 1], lw=1.5, color=c, solid_capstyle="round", zorder=8,
                path_effects=[pe.withStroke(linewidth=2.7, foreground="white")])
        for i in (len(curve) // 2, len(curve) - 6):
            if i + 4 >= len(curve): i = len(curve) - 8
            ax.annotate("", xy=curve[i + 4], xytext=curve[i],
                        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.24", color=c, lw=0,
                                        shrinkA=0, shrinkB=0, mutation_scale=11), zorder=9)
    for key, lab, col, off in [("GH", "Granulocyte lineage", "#2E6DA4", (0.4, 4.3)),
                               ("MH", "Monocytic/phagocytic lineage", "#7F4A12", (-5.4, -5.9)),
                               ("PH", "Prohemocyte", "#B2182B", (1.4, 4.6))]:
        m = (ct == "GH") if key == "GH" else ((ct == "MH2") | (ct == "TS")) if key == "MH" else (ct == "PH")
        mu = U[m].mean(0)
        ax.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.9, color=col, ha="center", va="center",
                zorder=10, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    if legend:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False, markerscale=4.5,
                  handletextpad=0.25, columnspacing=0.7, fontsize=6.0, labelspacing=0.25)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(xlabel + " 1"); ax.set_ylabel(xlabel + " 2")
    ax.set_box_aspect(1); ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
    ax.set_title(title, fontsize=7.5)
    p = os.path.join(FIG, out)
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

def main():
    dotplot(); cellcycle()
    G = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = G.obs["ct"].astype(str).replace({"MH1": "TS"}).values
    lineage(np.asarray(G.obsm["X_umap"]), ct, "Fig4E_lineages_UMAP_TS_preview.png", "UMAP", "Trajectory (UMAP)")
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    m = T.obs["vendor_cluster"].astype(str).map(GROUPS).isin(ORD).values
    ct2 = T.obs["vendor_cluster"].astype(str).map(GROUPS)[m].values
    lineage(np.asarray(T.obsm["X_tsne"])[m], ct2, "Fig4E_lineages_TSNE_TS_preview.png", "t-SNE", "Trajectory (t-SNE)")

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\blood_figs_ts.log", "w", encoding="utf8", buffering=1)
    main()
