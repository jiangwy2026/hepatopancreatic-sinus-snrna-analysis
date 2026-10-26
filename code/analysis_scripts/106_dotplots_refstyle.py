# -*- coding: utf-8 -*-
"""
106_dotplots_refstyle.py
(A) marker dot plot in the reference layout: genes on x (grouped by lineage with headers),
    subpopulations on y, size = pct.exp, colour = avg.exp.scale (z-score across subpopulations)
(B) cell-cycle dot plot: x = G1/S, G2, M ; y = subpopulations
"""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
pd.set_option("display.width", 250); sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ROWS = ["PH", "GH", "MH1", "MH2"]
CC = {"G1/S": ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"],
      "G2":   ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1"],
      "M":    ["KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]}

def dotpanel(ax, genes, rows, pct, z, smax=62, vmin=-2, vmax=2):
    ng = len(genes)
    for i in range(ng):
        for j in range(len(rows)):
            if np.isnan(pct[i, j]) or pct[i, j] <= 0.5:
                continue
            ax.scatter(i, j, s=3 + (pct[i, j] / 100) ** 1.25 * smax, c=[plt.cm.RdBu_r((np.clip(z[i, j], vmin, vmax) - vmin) / (vmax - vmin))],
                       edgecolors="#555555", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.7, ng - 0.3); ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xticks(range(ng)); ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=5.2, style="italic")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=7)
    for t in ax.get_yticklabels(): t.set_fontweight("bold")
    ax.tick_params(axis="both", length=0, pad=2)
    ax.grid(True, which="major", ls=(0, (1, 2)), lw=0.35, color="#CCCCCC", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
    ct = A.obs["vendor_cluster"].astype(str).map(base)
    H = A[ct.isin(ROWS).values].copy()
    H.obs["ct"] = pd.Categorical(ct[ct.isin(ROWS)].values, categories=ROWS)
    H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    sym2i = {}
    for i, s in enumerate(H.var["sym"].astype(str).values):
        sym2i.setdefault(s, i)
    print("cells:", H.shape, flush=True)

    def pctz(gl, sub=None):
        genes = [g for g in gl if g in sym2i]
        gi = [sym2i[g] for g in genes]
        X = H[:, gi].X
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        pct = np.zeros((len(genes), len(ROWS))); mean = np.zeros_like(pct)
        for j, g in enumerate(ROWS):
            m = (H.obs["ct"] == g).values
            Xs = X[m][:, :] if sub is None else None
            pct[:, j] = 100 * (X[m] > 0).mean(0)
            mean[:, j] = X[m].mean(0)
        z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)
        return genes, pct, z

    # ---------- A: marker dot plot (grouped by lineage)
    P = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"))
    blocks = [("Prohemocyte", P[P.grp == "PH"].sym.tolist()),
              ("Granulocyte", P[P.grp == "GH"].sym.tolist()),
              ("Monocytic hemocyte", P[P.grp.isin(["MH1", "MH2"])].sym.tolist())]
    genes, bounds = [], []
    for name, gl in blocks:
        genes += gl; bounds.append((name, len(genes) - len(gl), len(genes) - 1))
    g_ok, pct, z = pctz(genes)
    if len(g_ok) != len(genes):
        keep = [i for i, g in enumerate(genes) if g in g_ok]
        genes = [genes[i] for i in keep]; pct = pct[keep]; z = z[keep]
        for k, (name, a, b) in enumerate(bounds):
            bounds[k] = (name, sum(1 for i in keep if i < a), sum(1 for i in keep if i <= b) - 1)
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    dotpanel(ax, genes, ROWS, pct, z)
    for name, a, b in bounds:
        ax.plot([a - 0.45, b + 0.45], [-0.72, -0.72], color="black", lw=0.9, clip_on=False, zorder=4)
        ax.text((a + b) / 2, -1.05, name, ha="center", va="bottom", fontsize=6.4, weight="bold")
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + 0.42 * pos.height, 0.010, 0.34 * pos.height])
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-2, 2))
    cb = fig.colorbar(sm, cax=cax); cb.set_label("avg.exp.scale", fontsize=5.6, labelpad=1)
    cb.ax.tick_params(labelsize=5, length=1.2); cb.outline.set_visible(False)
    sx = pos.x1 + 0.055
    for k, p in enumerate((20, 40, 60, 80)):
        fig.add_artist(plt.Line2D([], []))
        ax.scatter([], [])
    for k, p in enumerate((20, 40, 60, 80)):
        y = 0.95 - k * 0.16
        fig.text(sx, pos.y0 + y * pos.height, "●", fontsize=4 + (p / 100) ** 1.25 * 62 ** 0.5, va="center", color="#666666")
        fig.text(sx + 0.030, pos.y0 + y * pos.height, "%d" % p, fontsize=5.2, va="center")
    fig.text(sx + 0.010, pos.y0 + 1.02 * pos.height, "pct.exp", fontsize=5.6, va="bottom")
    p1 = os.path.join(FIG, "Fig_blood4_dotplot_refstyle_preview.png")
    fig.savefig(p1, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p1, flush=True)

    # ---------- B: cell-cycle dot plot
    phase_genes = {k: [g for g in v if g in sym2i] for k, v in CC.items()}
    print("cell-cycle genes found:", {k: len(v) for k, v in phase_genes.items()}, flush=True)
    print({k: v for k, v in phase_genes.items()}, flush=True)
    ph_pct = np.zeros((len(CC), len(ROWS))); ph_z = np.zeros_like(ph_pct)
    for i, (ph, gl) in enumerate(phase_genes.items()):
        _, p, zz = pctz(gl)
        ph_pct[i] = p.mean(0); ph_z[i] = zz.mean(0)
    fig, ax = plt.subplots(figsize=(1.75, 2.5))
    dotpanel(ax, list(CC.keys()), ROWS, ph_pct, ph_z, smax=90)
    ax.set_title("Cell cycle", fontsize=7.5, pad=8)
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.06, pos.y0 + 0.55 * pos.height, 0.045, 0.22 * pos.height])
    cb = fig.colorbar(sm, cax=cax); cb.set_label("avg.exp.scale", fontsize=5.2, labelpad=1)
    cb.ax.tick_params(labelsize=4.6, length=1.1); cb.outline.set_visible(False)
    sx = pos.x1 + 0.10
    fig.text(sx, pos.y0 + 0.40 * pos.height, "pct.exp", fontsize=5.4)
    for k, p in enumerate((20, 40, 60, 80)):
        y = 0.28 - k * 0.085
        fig.text(sx, pos.y0 + y * pos.height, "●", fontsize=4 + (p / 100) ** 1.25 * 90 ** 0.5, va="center", color="#666666")
        fig.text(sx + 0.05, pos.y0 + y * pos.height, "%d" % p, fontsize=4.8, va="center")
    p2 = os.path.join(FIG, "Fig_cellcycle_dotplot_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p2, flush=True)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\dotplots.log", "w", encoding="utf8", buffering=1)
    main()
