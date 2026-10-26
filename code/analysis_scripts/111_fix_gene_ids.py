# -*- coding: utf-8 -*-
"""
111_fix_gene_ids.py - rebuild the 4-group dot plot and the single-gene t-SNE panels using the
gene IDs recorded by rank_genes_groups (symbols are ambiguous: e.g. 'Pxt' matches 8 gene IDs).
"""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
ROWS = ["PH", "GH", "MH1", "MH2"]

def dotpanel(ax, labels, rows, pct, z, smax=62, vmin=-2, vmax=2):
    ng = len(labels)
    for i in range(ng):
        for j in range(len(rows)):
            if pct[i, j] <= 0.5: continue
            ax.scatter(i, j, s=3 + (pct[i, j] / 100) ** 1.25 * smax,
                       c=[plt.cm.RdBu_r((np.clip(z[i, j], vmin, vmax) - vmin) / (vmax - vmin))],
                       edgecolors="#555555", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.7, ng - 0.3); ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xticks(range(ng)); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.2, style="italic")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=7)
    for t in ax.get_yticklabels(): t.set_fontweight("bold")
    ax.tick_params(axis="both", length=0, pad=2)
    ax.grid(True, which="major", ls=(0, (1, 2)), lw=0.35, color="#CCCCCC", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    ct_all = A.obs["vendor_cluster"].astype(str).map(BASE)
    H = A[ct_all.isin(ROWS).values].copy()
    H.obs["ct"] = pd.Categorical(ct_all[ct_all.isin(ROWS)].values, categories=ROWS)
    ids = list(H.var_names.astype(str)); sym = H.var["sym"].astype(str).values
    gpos = {g: i for i, g in enumerate(ids)}
    H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)

    P = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"))
    P["id"] = P["gid"].astype(str)
    gi = [gpos[g] for g in P.id]
    X = H[:, gi].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    pct = np.zeros((len(P), len(ROWS))); mean = np.zeros_like(pct)
    for j, g in enumerate(ROWS):
        m = (H.obs["ct"] == g).values
        pct[:, j] = 100 * (X[m] > 0).mean(0); mean[:, j] = X[m].mean(0)
    z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)
    out = pd.DataFrame(pct, index=P.sym.values, columns=ROWS).round(1)
    print("=== 修正后的 pct.exp（按基因 ID 取，行 = marker）===")
    print(out.to_string())
    out.to_csv(os.path.join(OUT, "blood4_dotplot_pct_fixed.csv"))

    # ---- corrected dot plot
    blocks = [("Prohemocyte", P[P.grp == "PH"]), ("Granulocyte", P[P.grp == "GH"]),
              ("Monocytic hemocyte", P[P.grp.isin(["MH1", "MH2"])])]
    order = np.concatenate([[k for k, r in enumerate(P.itertuples()) if r.grp in (blk[1].grp.values)] for blk in blocks])
    order = np.array([i for blk in blocks for i in blk[1].index])
    labels = P.sym.values[order]; pcto = pct[order]; zo = z[order]
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    dotpanel(ax, labels, ROWS, pcto, zo)
    k = 0
    for name, blk in blocks:
        a, b = k, k + len(blk) - 1; k = b + 1
        ax.plot([a - 0.45, b + 0.45], [-0.72, -0.72], color="black", lw=0.9, clip_on=False, zorder=4)
        ax.text((a + b) / 2, -1.05, name, ha="center", va="bottom", fontsize=6.4, weight="bold")
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
    p1 = os.path.join(FIG, "Fig_blood4_dotplot_refstyle_preview.png")
    fig.savefig(p1, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p1, flush=True)

    # ---- corrected single-gene t-SNE panels (all cells)
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    B = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))[T.obs_names].copy()
    ct = B.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    B.X = B.layers["counts"].copy(); sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
    bid = {g: i for i, g in enumerate(B.var_names.astype(str))}
    pick = ["Grik2", "ccnd2", "Pxt", "A2M", "PEN-3", "THBS4", "ITGA2B", "SDC2"]
    sel = P[P.sym.isin(pick)].drop_duplicates("sym").set_index("sym")
    genes = [(g, sel.loc[g, "grp"], sel.loc[g, "id"]) for g in pick if g in sel.index]
    print("\n=== 单基因面板用到的基因 ID ===")
    for g, grp, gid in genes: print("   %-7s %-4s %s" % (g, grp, gid), flush=True)
    rows = []
    for g, grp, gid in genes:
        v = np.asarray(B[:, bid[gid]].X.todense()).ravel()
        rows.append({"gene": g, "id": gid, "group": grp, "pct_all": round(100 * (v > 0).mean(), 1),
                     **{gg: round(100 * (v[ct == gg] > 0).mean(), 1) for gg in ROWS}})
    print(pd.DataFrame(rows).to_string(index=False))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "tsne_feature_gene_pct.csv"), index=False)

    fig, axes = plt.subplots(2, 4, figsize=(9.6, 5.1))
    for k, (g, grp, gid) in enumerate(genes):
        ax = axes.ravel()[k]
        v = np.asarray(B[:, bid[gid]].X.todense()).ravel()
        ax.scatter(tsne[:, 0], tsne[:, 1], c="#E8E8E8", s=0.8, linewidths=0, rasterized=True)
        pos = np.where(v > 0)[0]
        if len(pos):
            vmax = np.percentile(v[pos], 95); o = pos[np.argsort(-v[pos])]
            sc0 = ax.scatter(tsne[o, 0], tsne[o, 1], c=v[o], cmap="Reds", vmin=0, vmax=vmax,
                             s=1.8, linewidths=0, rasterized=True)
            cax = ax.inset_axes([0.87, 0.03, 0.10, 0.26])
            cb = fig.colorbar(sc0, cax=cax); cb.ax.tick_params(labelsize=4.2, length=1.0); cb.outline.set_visible(False)
        for gg in ROWS:
            m = ct == gg
            ax.text(np.median(tsne[m, 0]), np.median(tsne[m, 1]), gg, fontsize=5.4, color="#3A3A3A",
                    ha="center", va="center", zorder=5, path_effects=[pe.withStroke(linewidth=1.4, foreground="white")])
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_title("%s  (%s)   %0.1f%% cells" % (g, grp, 100 * (v > 0).mean()), fontsize=7.0, pad=3)
        for s_ in ax.spines.values(): s_.set_color("#999999")
    fig.suptitle("Single-gene expression across all %d cells (t-SNE; grey = not detected, dark red = high)" % T.n_obs,
                 fontsize=8, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    p2 = os.path.join(FIG, "Fig_tsne_genefeatures_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p2, flush=True)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\fix_ids.log", "w", encoding="utf8", buffering=1)
    main()
