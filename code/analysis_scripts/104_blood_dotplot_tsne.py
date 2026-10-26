# -*- coding: utf-8 -*-
"""
104_blood_dotplot_tsne.py
(1) recompute significant markers for the 4 blood groups (PH/GH/MH1/MH2) one-vs-rest
(2) dot plot of those markers
(3) t-SNE of the four groups
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
pd.set_option("display.width", 250); sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
MT = {"ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"}

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
    ct = A.obs["vendor_cluster"].astype(str).map(base)
    H = A[ct.isin(ORD).values].copy()
    H.obs["ct"] = pd.Categorical(ct[ct.isin(ORD)].values, categories=ORD)
    H.X = H.layers["counts"].copy()
    sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    symv = H.var["sym"].astype(str).values
    H.var["is_mt"] = np.isin(symv, list(MT))
    print("cells:", H.shape, "| groups:", H.obs["ct"].value_counts().to_dict(), flush=True)

    sc.tl.rank_genes_groups(H, "ct", method="wilcoxon", pts=True)
    mk = pd.concat([sc.get.rank_genes_groups_df(H, group=g).assign(grp=g) for g in ORD], ignore_index=True)
    mk["sym"] = mk["names"].map(dict(zip(H.var_names, symv))).fillna(mk["names"]).astype(str)
    mk["is_mt"] = mk["names"].map(H.var["is_mt"].to_dict()).fillna(False)
    mk["named"] = ~mk["sym"].str.startswith(("LOC", "gene-LOC"))
    mk.to_csv(os.path.join(OUT, "blood4_markers_all.csv"), index=False)
    import re
    HK = re.compile(r"^(Rp[LS]\d|rpl-|rps-|RpLP|eef|eIF|Act\d|Act5C|Tub|Hsp\d|hsc\d|HSC70|PRSS)", re.I)
    mk["hk"] = mk["sym"].str.match(HK).fillna(False)
    sig = mk[(~mk.is_mt) & (~mk.hk) & (mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)]
    print("genes significant (padj<0.05 & log2FC>0.5):")
    print(sig.groupby("grp", observed=True).size().to_string(), flush=True)

    TOP = int(os.environ.get("TOP", 8))
    picked, seen = [], set()
    for g in ORD:
        d = sig[sig.grp == g].sort_values(["named", "scores"], ascending=[False, False])
        cnt = 0
        for _, r in d.iterrows():
            if cnt >= TOP: break
            if r.sym in seen: continue
            picked.append({"grp": g, "sym": r.sym, "gid": r["names"], "log2FC": r.logfoldchanges, "q": r.pvals_adj})
            seen.add(r.sym); cnt += 1
    P = pd.DataFrame(picked)
    print("\nselected for display: %d genes" % len(P))
    print(P.groupby("grp", observed=True).sym.apply(lambda s: ", ".join(s)).to_string(), flush=True)
    P.to_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"), index=False)

    # ---- expression matrix for the dot plot
    idx = {g: i for i, g in enumerate(H.var_names.astype(str))}
    gi = [idx[g] for g in P.gid]
    X = H[:, gi].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    means = np.zeros((len(P), len(ORD))); pcts = np.zeros((len(P), len(ORD)))
    for j, g in enumerate(ORD):
        m = (H.obs["ct"] == g).values
        means[:, j] = X[m].mean(0)
        pcts[:, j] = 100 * (X[m] > 0).mean(0)
    Z = (means - means.min(1, keepdims=True)) / (means.max(1, keepdims=True) - means.min(1, keepdims=True) + 1e-9)

    # ---- dot plot
    fig, ax = plt.subplots(figsize=(2.4, 5.0))
    nG = len(P)
    for i in range(nG):
        for j in range(len(ORD)):
            ax.scatter(j, nG - 1 - i, s=4 + (pcts[i, j] / 100) ** 1.1 * 120, c=[plt.cm.magma(Z[i, j])],
                       edgecolors="#4D4D4D", linewidths=0.25, zorder=3)
    bounds = np.cumsum(P.groupby("grp", observed=True).size().reindex(ORD).fillna(0).values)
    for b in bounds[:-1]:
        ax.axhline(nG - b - 0.5, color="#BBBBBB", lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(-0.6, len(ORD) - 0.4); ax.set_ylim(-0.6, nG - 0.4)
    ax.set_yticks(range(nG)); ax.set_yticklabels([P.sym.iloc[nG - 1 - i] for i in range(nG)], fontsize=6.2)
    ax.set_xticks(range(len(ORD))); ax.set_xticklabels(ORD, fontsize=7)
    for t, g in zip(ax.get_xticklabels(), ORD): t.set_color(COL[g]); t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=1.5); ax.tick_params(axis="x", length=2, pad=1.5)
    for s in ("left", "bottom"): ax.spines[s].set_visible(False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_xlabel("Subpopulation", fontsize=7)
    # size legend + colourbar
    for k, pct in enumerate((25, 50, 75)):
        ax.scatter(len(ORD) - 0.25, nG - 1.2 - k * 1.35, s=4 + (pct / 100) ** 1.1 * 120, c="none",
                   edgecolors="#4D4D4D", linewidths=0.3)
        ax.text(len(ORD) - 0.05, nG - 1.2 - k * 1.35, "%d%%" % pct, fontsize=5.2, va="center")
    ax.text(len(ORD) - 0.25, nG + 0.1, "Expressing", fontsize=5.4, ha="center", va="bottom")
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.05, pos.y0 + 0.25 * pos.height, 0.045, 0.22 * pos.height])
    sm = plt.cm.ScalarMappable(cmap="magma", norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax); cb.set_label("Relative expression", fontsize=5.6)
    cb.ax.tick_params(labelsize=5, length=1.2); cb.outline.set_visible(False)
    p1 = os.path.join(FIG, "Fig_blood4_dotplot_preview.png")
    fig.savefig(p1, dpi=400, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
    print("saved", p1, flush=True)

    # ---- t-SNE of the four groups
    T = H.copy()
    sc.pp.highly_variable_genes(T, n_top_genes=2000)
    T = T[:, T.var.highly_variable].copy()
    sc.pp.scale(T, max_value=10)
    sc.tl.pca(T, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(T, n_neighbors=30, use_rep="X_pca")
    sc.tl.tsne(T, n_pcs=30, perplexity=30, random_state=0, n_jobs=8)
    H.obsm["X_tsne_blood"] = T.obsm["X_tsne"]
    H.write_h5ad(os.path.join(OUT, "adata_blood4_tsne.h5ad"))
    U = np.asarray(H.obsm["X_tsne_blood"])
    fig, ax = plt.subplots(figsize=(3.6, 3.5))
    for g in ORD:
        m = (H.obs["ct"] == g).values
        ax.scatter(U[m, 0], U[m, 1], s=1.4, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())))
    for g in ORD:
        m = (H.obs["ct"] == g).values
        cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
        ax.text(cx, cy, g, fontsize=7.5, weight="bold", ha="center", va="center", color="black",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=5,
              handletextpad=0.25, columnspacing=0.7, fontsize=6, labelspacing=0.25)
    p2 = os.path.join(FIG, "Fig_blood4_tsne_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
    print("saved", p2, flush=True)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\blood4_fig.log", "w", encoding="utf8", buffering=1)
    main()
