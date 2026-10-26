# -*- coding: utf-8 -*-
"""
139_annotation_evidence.py
注释依据图：(A) UMAP + 每个群的最佳匹配发表亚型与 z 值
(B) eLife 血细胞 6 亚型 × 本数据 10 群的 z-score 热图
(C) GEO 肝胰腺图谱 × 本数据 10 群的 z-score 热图
"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import TwoSlopeNorm
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
CL = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
NAME = {"0": "MH1", "1": "GH", "2": "PH", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NEU", "9": "NEU"}
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    U = np.asarray(A.obsm["X_umap"])
    cl = A.obs["vendor_cluster"].astype(str).values
    E = pd.read_csv(os.path.join(OUT, "elife_label_transfer_z.csv"), index_col=0)
    G = pd.read_csv(os.path.join(OUT, "geo_label_transfer_z.csv"), index_col=0)
    E.columns = [str(c) for c in E.columns]; G.columns = [str(c) for c in G.columns]
    EN = {"PH1(原血细胞1)": "PH1 (prohemocyte 1)", "PH2(原血细胞2)": "PH2 (prohemocyte 2)",
          "GH1(颗粒1)": "GH1 (granular 1)", "GH2(颗粒2)": "GH2 (granular 2)",
          "MH1(单核1)": "MH1 (monocytic 1)", "MH2(单核2·巨噬样)": "MH2 (monocytic 2, macrophage-like)"}
    E.index = [EN.get(str(i), str(i)) for i in E.index]
    XLAB = ["%s (C%s)" % (NAME[c], c) for c in CL]
    print("eLife 亚型:", list(E.index)); print("GEO 群:", list(G.index)[:6], "...", len(G.index))

    # 每群最佳匹配
    best = {}
    for c in CL:
        e = E[c] if c in E.columns else None; g = G[c] if c in G.columns else None
        eb = (e.idxmax(), e.max()) if e is not None else ("-", np.nan)
        gb = (g.idxmax(), g.max()) if g is not None else ("-", np.nan)
        best[c] = (eb, gb)
        print("C%s(%s): eLife %s %+.2f | GEO %s %+.2f" % (c, NAME[c], eb[0], eb[1], gb[0], gb[1]))

    fig = plt.figure(figsize=(11.6, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.15], wspace=0.32, left=0.03, right=0.985, top=0.88, bottom=0.16)
    # ---------- A: UMAP + 标注
    ax = fig.add_subplot(gs[0])
    for c in CL:
        m = cl == c
        ax.scatter(U[m, 0], U[m, 1], s=1.2, c=COL[NAME[c]], linewidths=0, rasterized=True)
    seen = set()
    for c in CL:
        m = cl == c
        cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
        g = NAME[c]
        if g in seen: continue
        seen.add(g)
        eb, gb = best[c]
        lab = g
        if eb[1] > 1.0: lab += "\n→ %s (%0.1f)" % (str(eb[0]).split(" (")[0], eb[1])
        elif gb[1] > 1.0: lab += "\n→ %s (%0.1f)" % (gb[0], gb[1])
        ax.text(cx, cy, lab, fontsize=5.6, ha="center", va="center", color="black", zorder=6,
                path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_box_aspect(1)
    ax.set_title("Our clusters → best published match (z)", fontsize=8)
    for s in ax.spines.values(): s.set_color("#999999")

    # ---------- B/C: z-score 热图
    def heat(ax, df, title, cols, xlabels):
        M = df[cols].astype(float).values
        vmax = np.nanmax(np.abs(M))
        im = ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax), aspect="auto")
        ax.set_xticks(range(len(xlabels))); ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=6)
        ax.set_yticks(range(len(df))); ax.set_yticklabels(df.index, fontsize=5.6)
        ax.set_title(title, fontsize=8)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if M[i, j] > 1.5:
                    ax.text(j, i, "%.1f" % M[i, j], fontsize=4.6, ha="center", va="center", color="white")
        return im
    axb = fig.add_subplot(gs[1])
    im1 = heat(axb, E, "Published hemocyte atlas\n(Yang et al. 2022, eLife)", CL, XLAB)
    axc = fig.add_subplot(gs[2])
    im2 = heat(axc, G, "Published hepatopancreas atlas\n(GSE306674)", CL, XLAB)
    cax = fig.add_axes([0.955, 0.25, 0.008, 0.55])
    cb = fig.colorbar(im2, cax=cax); cb.set_label("marker-set enrichment (z)", fontsize=6)
    cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
    fig.suptitle("Annotation evidence: label transfer of published marker sets (columns = our clusters)", fontsize=8.5, y=0.99)
    p = os.path.join(FIG, "Fig_annotation_evidence_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\annotation_evidence.log", "w", encoding="utf8", buffering=1)
    main()
