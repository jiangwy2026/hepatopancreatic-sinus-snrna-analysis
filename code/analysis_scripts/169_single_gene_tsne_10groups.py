# -*- coding: utf-8 -*-
"""169_single_gene_tsne_10groups.py - 每组 4 个 marker，共 40 格单基因 t-SNE"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
GROUPS = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER = ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
COL = {"PH":"#B2182B","GH":"#E08214","TS":"#D6604D","MH2":"#8C510A","F":"#006D2C","R":"#41AB5D",
       "B":"#78C679","M":"#2166AC","MyoEpi":"#92C5DE","PM-like":"#7B3294"}
GREY = "#C9C9C9"
GCMAP = ListedColormap(plt.cm.Greens(np.linspace(0.28, 0.95, 256)))

def main():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))[T.obs_names].copy()
    ct = A.obs["vendor_cluster"].astype(str).map(GROUPS).values
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}

    M = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv"))
    M["named"] = ~M["sym"].astype(str).str.match(r"^(LOC|gene-LOC)")
    picks = []
    for g in ORDER:
        d = M[(M.grp == g) & (M["rank"] <= 10)]
        sel = list(d.sort_values("rank").head(2).itertuples())
        for r in d[d.named].itertuples():
            if len(sel) >= 4: break
            if r.sym not in [x.sym for x in sel]: sel.append(r)
        for r in d.itertuples():
            if len(sel) >= 4: break
            if r.sym not in [x.sym for x in sel]: sel.append(r)
        picks += [(g, x.sym, x.gid) for x in sel[:4]]
    print("每群选基因:", {g: [p[1] for p in picks if p[0] == g] for g in ORDER})

    ncol, nrow = 8, 5
    fig = plt.figure(figsize=(15.0, 9.4))
    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.46, wspace=0.40,
                  left=0.02, right=0.975, top=0.93, bottom=0.03)
    for k, (g, sym, gid) in enumerate(picks):
        r, c = divmod(k, ncol)
        ax = fig.add_subplot(gs[r, c])
        v = np.asarray(A[:, gpos[gid]].X.todense()).ravel()
        ax.scatter(tsne[:, 0], tsne[:, 1], c=GREY, s=0.7, linewidths=0, rasterized=True)
        pos = np.where(v > 0)[0]
        if len(pos):
            vmax = np.percentile(v[pos], 95)
            rel = np.clip(v[pos] / vmax, 0, 1); o = np.argsort(rel)
            sm = ax.scatter(tsne[pos[o], 0], tsne[pos[o], 1], c=rel[o], cmap=GCMAP, vmin=0, vmax=1,
                            s=1.1, linewidths=0, rasterized=True)
            cax = ax.inset_axes([1.035, 0.30, 0.045, 0.40])
            cb = fig.colorbar(sm, cax=cax); cb.set_ticks([0, 1]); cb.ax.tick_params(labelsize=3.4, length=0.7, pad=0.8)
            cb.outline.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        in_g = v[ct == g]
        ax.set_title("%s [%s]\n%.0f%% of group / %.1f%% all" % (sym, g, 100 * (in_g > 0).mean(), 100 * (v > 0).mean()),
                     fontsize=4.7, pad=2)
        for s_ in ax.spines.values(): s_.set_color("#BBBBBB"); s_.set_linewidth(0.4)
    fig.suptitle("Single-gene expression across all %d cells (t-SNE; grey = not detected, dark green = high)"
                 % T.n_obs, fontsize=9, y=0.985)
    p = os.path.join(FIG, "Fig_tsne_genefeatures_10groups_preview.png")
    fig.savefig(p, dpi=350, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)
    pd.DataFrame(picks, columns=["group", "gene", "gid"]).to_csv(os.path.join(OUT, "single_gene_panel_10groups.csv"), index=False)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\single_gene_10.log", "w", encoding="utf8", buffering=1)
    main()
