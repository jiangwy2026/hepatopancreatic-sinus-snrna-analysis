# -*- coding: utf-8 -*-
"""167_fig1_fig2_10groups.py - 按最终 10 群重出 Fig1(UMAP)、Fig1b(t-SNE)、Fig2(marker 热图)"""
import io, sys, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
GROUPS = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER = ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
FULL = {"PH":"Prohemocyte","GH":"Granular hemocyte","TS":"Biosynthetic transition state","MH2":"Macrophage-like phagocyte",
        "F":"F cell","R":"R cell","B":"B cell","M":"M cell","MyoEpi":"Myoepithelial-like","PM-like":"Peritrophic matrix-like"}
COL = {"PH":"#B2182B","GH":"#E08214","TS":"#D6604D","MH2":"#8C510A","F":"#006D2C","R":"#41AB5D",
       "B":"#78C679","M":"#2166AC","MyoEpi":"#92C5DE","PM-like":"#7B3294"}

def label_panel(ax, U, ct, order, fs=7.5):
    for g in order:
        m = ct == g
        if m.sum() == 0: continue
        cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
        ax.text(cx, cy, g, fontsize=fs, weight="bold", ha="center", va="center", color="black",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")], zorder=6)

def main():
    H = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    ct = H.obs["vendor_cluster"].astype(str).map(GROUPS).values
    U = np.asarray(H.obsm["X_umap"]); T = np.asarray(H.obsm["X_tsne"])

    # ---- Fig1 UMAP
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for g in ORDER:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=2.0, c=COL[g], linewidths=0, rasterized=True,
                   label="%s  %s (%d)" % (g, FULL[g], int(m.sum())))
    label_panel(ax, U, ct, ORDER)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_box_aspect(1)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
              fontsize=5.6, handletextpad=0.3, labelspacing=0.3)
    p = os.path.join(FIG, "Fig1_UMAP_10groups_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # ---- Fig1b t-SNE
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for g in ORDER:
        m = ct == g
        ax.scatter(T[m, 0], T[m, 1], s=2.0, c=COL[g], linewidths=0, rasterized=True,
                   label="%s  %s (%d)" % (g, FULL[g], int(m.sum())))
    label_panel(ax, T, ct, ORDER)
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
              fontsize=5.6, handletextpad=0.3, labelspacing=0.3)
    p = os.path.join(FIG, "Fig1b_TSNE_10groups_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # ---- Fig2 heatmap（单细胞级，宽度按 log 细胞数）
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.obs["ct"] = pd.Categorical(A.obs["vendor_cluster"].astype(str).map(GROUPS), categories=ORDER)
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    P = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv"))
    P = P[P["rank"] <= 5].copy()
    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
    rng = np.random.default_rng(0); cell_idx, marks, shown, total = [], [], [], []
    for g in ORDER:
        idx = np.where((A.obs["ct"] == g).values)[0]
        n = len(idx); k = min(n, max(15, int(round(25 * math.log10(n)))))
        pick = np.sort(rng.choice(idx, k, replace=False))
        cell_idx += list(pick); marks.append(len(cell_idx)); shown.append(k); total.append(n)
    cell_idx = np.array(cell_idx)
    gi = [gpos[g] for g in P["gid"].astype(str)]
    E = A[cell_idx][:, gi].X.toarray().astype(np.float32)
    Z = np.clip(E / (np.percentile(E, 99, axis=0) + 1e-9), 0, 1).T
    fig = plt.figure(figsize=(9.4, 5.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.030, 0.020, 1.0], width_ratios=[0.26, 1.0, 0.026],
                          hspace=0.02, wspace=0.03, left=0.02, right=0.965, top=0.965, bottom=0.03)
    axb = fig.add_subplot(gs[0, 1]); axn = fig.add_subplot(gs[1, 1])
    axh = fig.add_subplot(gs[2, 1]); axl = fig.add_subplot(gs[2, 0]); axc = fig.add_subplot(gs[2, 2])
    nC = Z.shape[1]
    bars = np.zeros((1, nC))
    for i in range(len(ORDER)):
        a = 0 if i == 0 else marks[i - 1]; bars[0, a:marks[i]] = i
    axb.imshow(bars, aspect="auto", cmap=ListedColormap([COL[g] for g in ORDER]),
               norm=BoundaryNorm(np.arange(-0.5, len(ORDER) + 0.5), len(ORDER)), extent=(0, nC, 0, 1))
    axb.set_xticks([]); axb.set_yticks([])
    for i, g in enumerate(ORDER):
        a = 0 if i == 0 else marks[i - 1]; b = marks[i]
        axb.text((a + b) / 2, 0.5, g, ha="center", va="center", fontsize=5.0, weight="bold", color="white")
        if i > 0: axb.axvline(a, color="white", lw=0.8)
    for sp in axb.spines.values(): sp.set_visible(False)
    axn.set_xlim(0, nC); axn.set_ylim(0, 1); axn.axis("off")
    for i, g in enumerate(ORDER):
        a = 0 if i == 0 else marks[i - 1]; b = marks[i]
        axn.text((a + b) / 2, 0.0, "n=%s" % format(total[i], ","), ha="center", va="bottom", fontsize=3.8)
    im = axh.imshow(Z, aspect="auto", cmap="magma", vmin=0, vmax=1, interpolation="nearest")
    axh.set_xticks([]); axh.set_yticks([])
    for sp in axh.spines.values(): sp.set_visible(False)
    for i in range(1, len(ORDER)): axh.axvline(marks[i - 1], color="white", lw=0.8)
    axl.set_xlim(0, 1); axl.set_ylim(Z.shape[0] - 0.5, -0.5); axl.axis("off")
    for i, r in enumerate(P.itertuples()):
        axl.text(1.0, i, "%s (%s)" % (r.sym, r.gid), ha="right", va="center", fontsize=4.4)
    pos = axh.get_position(); cw, ch = 0.012, pos.height * 0.18
    axc.set_position([pos.x1 + 0.006, pos.y1 - ch - 0.012, cw, ch])
    cb = fig.colorbar(im, cax=axc); cb.set_label("Relative expression", fontsize=4.8)
    cb.ax.tick_params(labelsize=4.2, length=1.1); cb.outline.set_visible(False)
    p = os.path.join(FIG, "Fig2_marker_heatmap_10groups_preview.png")
    fig.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\fig12_10groups.log", "w", encoding="utf8", buffering=1)
    main()
