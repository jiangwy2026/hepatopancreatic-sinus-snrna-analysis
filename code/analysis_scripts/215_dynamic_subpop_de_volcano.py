# -*- coding: utf-8 -*-
"""215_dynamic_subpop_de_volcano.py
Within-subpopulation differential expression: Hemospasia vs Control.
Produces five volcano panels and DE tables for the dynamic subpopulations.
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
from matplotlib.lines import Line2D

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-20", "第二张图_抽血处理细胞变化")
DE_DIR = os.path.join(OUT, "dynamic_subpop_de")
os.makedirs(PARTS, exist_ok=True); os.makedirs(DE_DIR, exist_ok=True); os.makedirs(DELIV, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style

apply_pub_style(7)
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
sc.settings.verbosity = 0

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
DYN = [("PH", 2.30, "up"), ("MH2", 2.17, "up"), ("R", 2.29, "up"),
       ("MyoEpi", 0.44, "down"), ("PM-like", 0.27, "down")]
UP = "#B2182B"; DOWN = "#2166AC"; NS = "#BDBDBD"; INK = "#222222"


def run_de():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.obs["ct"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    sym = dict(zip(A.var_names.astype(str), A.var["sym"].astype(str)))
    mt = dict(zip(A.var_names.astype(str), A.var["mt"].astype(bool)))
    tables = {}
    for grp, fold, direction in DYN:
        sub = A[A.obs["ct"].astype(str) == grp].copy()
        print(f"DE {grp}: {sub.n_obs} cells; group counts {sub.obs['group'].astype(str).value_counts().to_dict()}")
        sc.tl.rank_genes_groups(sub, "group", groups=["treated"], reference="blank",
                                method="wilcoxon", pts=True, tie_correct=False)
        d = sc.get.rank_genes_groups_df(sub, group="treated")
        d["gid"] = d["names"].astype(str)
        d["sym"] = d["gid"].map(sym)
        d["is_mt"] = d["gid"].map(mt).fillna(False)
        d = d[~d["is_mt"]].copy()
        d = d[np.isfinite(d["logfoldchanges"]) & np.isfinite(d["pvals_adj"])].copy()
        d["neglogq"] = -np.log10(np.clip(d["pvals_adj"].astype(float), 1e-300, 1.0))
        d["direction"] = np.where(d["logfoldchanges"] > 0, "up", "down")
        d = d.sort_values(["pvals_adj", "scores"], ascending=[True, False])
        d.to_csv(os.path.join(DE_DIR, f"DE_{grp}_treated_vs_control.csv"), index=False, encoding="utf-8-sig")
        tables[grp] = d.copy()
    return tables


def pick_labels(d, grp):
    sig = d[(d.pvals_adj < 0.05) & (d.logfoldchanges.abs() > 0.5)].copy()
    if sig.empty:
        sig = d.head(10).copy()
    sig["abs_lfc"] = sig.logfoldchanges.abs()
    sig = sig.sort_values(["neglogq", "abs_lfc"], ascending=False)
    labels = []
    seen = set()
    for r in sig.itertuples():
        s = str(r.sym)
        if s not in seen:
            labels.append(r._asdict())
            seen.add(s)
        if len(labels) >= 3:
            break
    if grp == "PH" and "Grik2" in set(d.sym.astype(str)):
        g = d[d.sym.astype(str) == "Grik2"]
        if len(g) and "Grik2" not in seen:
            labels.append(g.iloc[0].to_dict())
    return labels


def main():
    tables = run_de()
    fig = plt.figure(figsize=(7.2, 6.6), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=0.075, right=0.975, top=0.90, bottom=0.150,
                          wspace=0.34, hspace=0.48)
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    for (grp, fold, direction), pos in zip(DYN, positions):
        ax = fig.add_subplot(gs[pos])
        d = tables[grp]
        sig = (d.pvals_adj < 0.05) & (d.logfoldchanges.abs() > 0.5)
        up = sig & (d.logfoldchanges > 0)
        dn = sig & (d.logfoldchanges < 0)
        ns = ~sig
        ax.scatter(d.loc[ns, "logfoldchanges"], np.minimum(d.loc[ns, "neglogq"], 100),
                   s=3, c=NS, alpha=0.35, linewidths=0, rasterized=True, zorder=2)
        ax.scatter(d.loc[dn, "logfoldchanges"], np.minimum(d.loc[dn, "neglogq"], 100),
                   s=4, c=DOWN, alpha=0.65, linewidths=0, rasterized=True, zorder=3)
        ax.scatter(d.loc[up, "logfoldchanges"], np.minimum(d.loc[up, "neglogq"], 100),
                   s=4, c=UP, alpha=0.65, linewidths=0, rasterized=True, zorder=3)
        ax.axhline(-np.log10(0.05), color="#777777", lw=0.45, ls=(0, (3, 2)))
        ax.axvline(-0.5, color="#777777", lw=0.45, ls=(0, (3, 2)))
        ax.axvline(0.5, color="#777777", lw=0.45, ls=(0, (3, 2)))
        ax.set_xlim(-4.5, 6.5)
        ax.set_ylim(0, 105)
        ax.set_xlabel("log$_2$ fold change (Hemospasia / Control)", fontsize=5.5)
        ax.set_ylabel("−log$_{10}$ adjusted $P$", fontsize=5.5)
        ax.tick_params(labelsize=5.0, length=1.5, width=0.45)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="both", color="#F0F0F0", lw=0.3, zorder=0)
        col = UP if direction == "up" else DOWN
        ax.set_title(f"{grp}  {'↑' if direction == 'up' else '↓'} {fold:.2f}×",
                     fontsize=6.2, color=col, fontweight="bold", pad=4)
        # label selected genes
        for k, r in enumerate(pick_labels(d, grp)):
            x = float(r["logfoldchanges"]); y = min(float(r["neglogq"]), 100)
            dx = 0.15 if x >= 0 else -0.15
            ax.annotate(str(r["sym"]), xy=(x, y), xytext=(x + dx, min(y + 3.0, 101)),
                        fontsize=5.0, color=INK, ha="left" if x >= 0 else "right",
                        va="bottom", arrowprops=dict(arrowstyle="-", lw=0.3, color="#777777"))
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=UP, markersize=3.2,
                      label="up in Hemospasia"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=DOWN, markersize=3.2,
                      label="down in Hemospasia"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=NS, markersize=3.2,
                      label="not passing thresholds")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.050), ncol=3,
               frameon=False, fontsize=5.2, handlelength=0.7, columnspacing=1.2)
    fig.text(0.50, 0.018,
             "Thresholds: adjusted P < 0.05 and |log2FC| > 0.5. Cell-level Wilcoxon; "
             "no biological replicates.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    fig.text(0.50, 0.965, "Within-subpopulation differential expression: Hemospasia vs Control",
             ha="center", va="top", fontsize=7.5)
    out = os.path.join(PARTS, "Fig5_dynamic_subpops_DE_volcano_preview.png")
    fig.savefig(out, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(os.path.join(DELIV, "FigS5_dynamic_DE_volcano.svg"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(os.path.join(DELIV, "FigS5_dynamic_DE_volcano.tiff"), dpi=1200, bbox_inches="tight", pad_inches=0.04, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    print("saved", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
