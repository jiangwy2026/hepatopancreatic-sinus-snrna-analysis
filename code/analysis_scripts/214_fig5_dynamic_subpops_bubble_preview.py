# -*- coding: utf-8 -*-
"""214_fig5_dynamic_subpops_bubble_preview.py
Preview bubble panels for composition-changing subpopulations.

Rows: selected marker genes; columns: Control and Hemospasia.
Dot size: percent expressing. Colour: mean log-normalized expression.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
os.makedirs(PARTS, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style

apply_pub_style(7)
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
DYN = {
    "PH": ["PPAF2", "Igfbp7", "SSPO", "Grik2", "ccnd2", "LOC113808957", "slo", "LOC113823945"],
    "MH2": ["LOC113810024", "TIE1", "ITGA2B", "PTPRA", "SDC2", "Tl", "AGAP007347", "LOC113811918"],
    "R": ["mesh", "Oacyl", "ACE", "HSPG2", "hdc", "LOC113802931", "LOC113823669", "LOC113815680"],
    "MyoEpi": ["LOC113814314", "lcc1", "A1i3", "KCNQ5", "MINPP1"],
    "PM-like": ["LOC113809538", "LOC113809539", "Fas2", "LOC113809540", "GPx", "hsc71", "LOC113829816", "LOC113805076"],
}
LAYOUT = [("PH", 0, 0), ("MH2", 0, 1), ("R", 0, 2), ("MyoEpi", 1, 0), ("PM-like", 1, 1)]
UP = "#B2182B"; DOWN = "#2166AC"


def resolve_gid(M, grp, gene):
    sub = M[(M.grp == grp) & (M.sym == gene)].sort_values("scores", ascending=False)
    if sub.empty:
        raise ValueError(f"Missing marker {grp} {gene}")
    return str(sub.iloc[0]["names"])


def main():
    M = pd.read_csv(os.path.join(OUT, "markers_10groups.csv"))
    comp = pd.read_csv(os.path.join(OUT, "composition_new_names_10groups.csv"))
    comp = comp.set_index("group")
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.obs["ct"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    A.obs["grp2"] = A.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    rows, gene_meta = [], []
    for grp in DYN:
        sub = A[A.obs["ct"].astype(str) == grp]
        for rank, gene in enumerate(DYN[grp], 1):
            gid = resolve_gid(M, grp, gene)
            if gid not in set(A.var_names.astype(str)):
                raise ValueError(f"Gene id not in AnnData: {gid}")
            x = np.asarray(sub[:, gid].X.todense()).ravel().astype(float)
            gene_meta.append({"group": grp, "rank": rank, "gene": gene, "gid": gid})
            for group in ("Control", "Hemospasia"):
                m = (sub.obs["grp2"].astype(str).values == group)
                rows.append({"group": grp, "gene": gene, "gid": gid,
                             "sample_group": group,
                             "pct_expressing": float(100 * (x[m] > 0).mean()),
                             "mean_logexpr": float(x[m].mean())})
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(OUT, "Fig5_dynamic_subpops_bubble_data.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(gene_meta).to_csv(os.path.join(OUT, "Fig5_dynamic_subpops_bubble_genes.csv"),
                                   index=False, encoding="utf-8-sig")

    vmax = max(1.0, float(D["mean_logexpr"].quantile(0.99)))
    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=0, vmax=vmax)
    fig = plt.figure(figsize=(7.20, 4.90), facecolor="white")
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.55], height_ratios=[1, 1],
                          left=0.065, right=0.985, top=0.855, bottom=0.100,
                          wspace=0.58, hspace=0.42)
    for grp, r, c in LAYOUT:
        ax = fig.add_subplot(gs[r, c])
        genes = DYN[grp]
        n = len(genes)
        sub = D.loc[D.group == grp].copy()
        n_ctrl = int(comp.loc[grp, "n_blank"])
        n_trt = int(comp.loc[grp, "n_treated"])
        fold = float(comp.loc[grp, "fold_treated_vs_blank"])
        direction = "↑" if fold >= 1 else "↓"
        tcol = UP if fold >= 1 else DOWN
        for i, gene in enumerate(genes):
            for j, sample_group in enumerate(("Control", "Hemospasia")):
                r0 = sub[(sub.gene == gene) & (sub.sample_group == sample_group)].iloc[0]
                ax.scatter(j, i, s=6 + 58 * float(r0.pct_expressing) / 100.0,
                           c=[cmap(norm(float(r0.mean_logexpr)))],
                           edgecolors="white", linewidths=0.35, zorder=3)
        for i in range(n):
            ax.axhline(i, color="#ECECEC", lw=0.35, zorder=0)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(n - 0.55, -1.30)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Control", "Hemospasia"], fontsize=5.0)
        ax.set_yticks(np.arange(n))
        ax.set_yticklabels(genes, fontsize=5.0)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.text(0.50, 1.12, f"{grp}  {direction} {fold:.2f}×", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=tcol)
        ax.text(0.50, 1.03, f"n = {n_ctrl} / {n_trt} cells", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=5.0, color="#555555")

    # shared colour bar
    cax = fig.add_axes([0.825, 0.52, 0.020, 0.30])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("Mean log expression", fontsize=5.2, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.5, width=0.4)
    cb.outline.set_linewidth(0.4)

    # size legend
    axl = fig.add_axes([0.805, 0.14, 0.16, 0.22])
    axl.set_xlim(0, 1); axl.set_ylim(0, 1); axl.axis("off")
    axl.text(0.50, 0.95, "% expressing", ha="center", va="top", fontsize=5.0, color="#333333")
    for p, y in zip([25, 50, 75], [0.68, 0.43, 0.18]):
        axl.scatter([0.26], [y], s=6 + 58 * p / 100.0, c="white",
                    edgecolors="#555555", linewidths=0.5)
        axl.text(0.52, y, f"{p}%", ha="left", va="center", fontsize=5.0, color="#333333")

    fig.text(0.50, 0.985, "Dynamic subpopulations: marker expression by group",
             ha="center", va="top", fontsize=7.0)
    fig.text(0.50, 0.018,
             "Arrows/fold change: cell proportion, Hemospasia vs Control. "
             "Dot size: percent expressing; colour: mean log-normalized expression.",
             ha="center", va="bottom", fontsize=5.0, color="#444444")
    out = os.path.join(PARTS, "Fig5_dynamic_subpops_bubble_preview.png")
    fig.savefig(out, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
