# -*- coding: utf-8 -*-
"""
107_cellcycle_fig.py - cell-cycle figure: (A) reference-style dot plot (G1/S, G2, M x subtype)
(B) phase composition per subtype.  Dot sizes in A are scaled to the panel maximum because the
absolute detection of cell-cycle transcripts is low in this shallow library.
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
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ROWS = ["PH", "GH", "MH1", "MH2"]
CC = {"G1/S": ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"],
      "G2":   ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1"],
      "M":    ["KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]}
PHCOL = {"G1": "#BDBDBD", "S": "#4393C3", "G2M": "#B2182B"}

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
    ct = A.obs["vendor_cluster"].astype(str).map(base)
    H = A[ct.isin(ROWS).values].copy()
    H.obs["ct"] = pd.Categorical(ct[ct.isin(ROWS)].values, categories=ROWS)
    H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    sym2i = {}; sym2id = {}
    for i, s in enumerate(H.var["sym"].astype(str).values):
        sym2i.setdefault(s, i); sym2id.setdefault(s, H.var_names[i])

    ph_pct = np.zeros((len(CC), len(ROWS))); ph_z = np.zeros_like(ph_pct)
    for i, (ph, gl) in enumerate(CC.items()):
        genes = [g for g in gl if g in sym2i]
        X = H[:, [sym2i[g] for g in genes]].X
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        pct = np.zeros((len(genes), len(ROWS))); mean = np.zeros_like(pct)
        for j, g in enumerate(ROWS):
            m = (H.obs["ct"] == g).values
            pct[:, j] = 100 * (X[m] > 0).mean(0); mean[:, j] = X[m].mean(0)
        z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)
        ph_pct[i] = pct.mean(0); ph_z[i] = z.mean(0)
    print("pct.exp matrix (rows=phase, cols=group):"); print(pd.DataFrame(ph_pct, index=list(CC), columns=ROWS).round(1).to_string())
    print("avg.exp.scale:"); print(pd.DataFrame(ph_z, index=list(CC), columns=ROWS).round(2).to_string())

    fig = plt.figure(figsize=(5.6, 2.5))
    ax1 = fig.add_axes([0.10, 0.20, 0.30, 0.62]); ax2 = fig.add_axes([0.53, 0.20, 0.31, 0.62])
    smax_pct = max(45.0, np.nanmax(ph_pct))
    for i in range(len(CC)):
        for j in range(len(ROWS)):
            s = 3 + (ph_pct[i, j] / smax_pct) ** 1.15 * 150
            ax1.scatter(i, j, s=s, c=[plt.cm.RdBu_r((np.clip(ph_z[i, j], -2, 2) + 2) / 4)],
                        edgecolors="#555555", linewidths=0.25, zorder=3)
    ax1.set_xlim(-0.7, len(CC) - 0.3); ax1.set_ylim(len(ROWS) - 0.45, -0.55)
    ax1.set_xticks(range(len(CC))); ax1.set_xticklabels(list(CC.keys()), fontsize=7)
    ax1.set_yticks(range(len(ROWS))); ax1.set_yticklabels(ROWS, fontsize=7.5, fontweight="bold")
    ax1.tick_params(length=0, pad=2); ax1.grid(True, ls=(0, (1, 2)), lw=0.35, color="#CCCCCC", zorder=0)
    for s_ in ax1.spines.values(): s_.set_visible(False)
    ax1.set_title("Cell cycle genes", fontsize=7.5, pad=6)
    for k, p in enumerate((10, 20, 30, 40)):
        ax1.scatter(2.55, -0.35 + k * 0.55, s=3 + (p / smax_pct) ** 1.15 * 150, c="none",
                    edgecolors="#666666", linewidths=0.3, clip_on=False)
        ax1.text(2.75, -0.35 + k * 0.55, "%d%%" % p, fontsize=5, va="center", clip_on=False)
    ax1.text(2.55, -0.75, "pct.exp", fontsize=5.4, ha="center", va="bottom", clip_on=False)
    fig.text(0.435, 0.86, "avg.exp.scale", fontsize=5.6, ha="center")
    cax = fig.add_axes([0.415, 0.42, 0.012, 0.30])
    from matplotlib.cm import ScalarMappable
    cb = fig.colorbar(ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-2, 2)), cax=cax)
    cb.ax.tick_params(labelsize=4.8, length=1.1); cb.outline.set_visible(False)

    comp = pd.read_csv(os.path.join(OUT, "cellcycle_phase_composition.csv"), index_col=0).reindex(ROWS)
    order = ["G1", "S", "G2M"]
    left = np.zeros(len(ROWS))
    for ph in order:
        v = comp[ph].values
        ax2.barh(np.arange(len(ROWS)), v, left=left, color=PHCOL[ph], edgecolor="white", height=0.62,
                 label=ph)
        for k, (l, x) in enumerate(zip(left, v)):
            if x >= 6:
                ax2.text(l + x / 2, k, "%.0f" % x, fontsize=5.4, ha="center", va="center",
                         color="white" if ph != "G1" else "#333333")
        left += v
    ax2.set_xlim(0, 100); ax2.set_ylim(len(ROWS) - 0.45, -0.55)
    ax2.set_yticks(range(len(ROWS))); ax2.set_yticklabels([])
    ax2.set_xlabel("Cells per phase (%)", fontsize=6.5)
    ax2.tick_params(axis="x", labelsize=5.5, length=1.6, width=0.5)
    for s_ in ("top", "right", "left"): ax2.spines[s_].set_visible(False)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=3, frameon=False, fontsize=5.8,
               handlelength=0.9, handletextpad=0.3, columnspacing=0.8)
    ax2.set_title("Cell-cycle phase composition", fontsize=7.5, pad=16)
    fig.text(0.055, 0.93, "A", fontsize=9, fontweight="bold")
    fig.text(0.485, 0.93, "B", fontsize=9, fontweight="bold")
    p = os.path.join(FIG, "Fig_cellcycle_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\cellcycle_fig.log", "w", encoding="utf8", buffering=1)
    main()
