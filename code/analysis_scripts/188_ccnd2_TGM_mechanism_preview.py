# -*- coding: utf-8 -*-
"""188_ccnd2_TGM_mechanism_preview.py
Early-PH mechanism figure: volcano, ccnd2/TGM expression, and candidate model.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import mannwhitneyu

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
CTRL = "#4C78A8"; TRT = "#D6604D"


def draw_volcano(ax):
    D = pd.read_csv(os.path.join(OUT, "PH_early_treated_vs_blank_DE.csv"))
    D = D[np.isfinite(D.logfoldchanges) & np.isfinite(D.pvals_adj)].copy()
    D["neglogq"] = -np.log10(np.clip(D.pvals_adj.values, 1e-300, 1))
    sig = (D.pvals_adj < 0.05) & (D.logfoldchanges.abs() > 0.5)
    up = sig & (D.logfoldchanges > 0); dn = sig & (D.logfoldchanges < 0)
    ax.scatter(D.loc[~sig, "logfoldchanges"], D.loc[~sig, "neglogq"], s=2.4, c="#BDBDBD",
               linewidths=0, alpha=0.55, rasterized=True)
    ax.scatter(D.loc[up, "logfoldchanges"], D.loc[up, "neglogq"], s=5.2, c="#D73027",
               linewidths=0, alpha=0.88, rasterized=True)
    ax.scatter(D.loc[dn, "logfoldchanges"], D.loc[dn, "neglogq"], s=5.2, c="#4575B4",
               linewidths=0, alpha=0.75, rasterized=True)
    ax.axvline(-0.5, color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax.axvline(0.5, color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax.axhline(-np.log10(0.05), color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax.set_xlim(-6.5, 3.5); ax.set_ylim(0, 30)
    ax.set_xlabel("log$_2$ fold change (treated / control)", fontsize=5.8)
    ax.set_ylabel("−log$_{10}$ adjusted $P$", fontsize=5.8)
    ax.tick_params(labelsize=5.0, length=1.2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.set_title("Early PH response", fontsize=7.0, pad=4)
    labels = {
        "LOC113823945": ("hemocyte TGM", (-1.3, 2.0), "#B2182B"),
        "ccnd2": ("ccnd2", (1.0, 1.0), "#B2182B"),
        "SSPO": ("SSPO", (-0.2, 1.4), "#B2182B"),
        "LOC113808957": ("Notch-like", (0.4, 2.8), "#B2182B"),
        "THBS4": ("THBS4", (-0.8, 3.2), "#2166AC"),
        "PEN-3": ("PEN-3", (-1.0, 2.6), "#2166AC"),
        "LOC113801825": ("crustin", (-0.8, 2.4), "#2166AC"),
        "LOC113810905": ("WAP", (-1.4, 2.4), "#2166AC"),
    }
    for gene, (lab, off, col) in labels.items():
        r = D[D.gene == gene]
        if r.empty: continue
        r = r.sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False]).iloc[0]
        y = min(r.neglogq, 29.5)
        ax.annotate(lab, xy=(r.logfoldchanges, y),
                    xytext=(r.logfoldchanges + off[0], min(y + off[1], 29.5)),
                    fontsize=5.1, color=col, ha="left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.35, color=col, alpha=0.8))


def draw_violin(ax):
    A = ad.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    A.obs["pt_bin3"] = pd.qcut(A.obs["dpt_pseudotime"].astype(float), 3, labels=["P1", "P2", "P3"])
    early = A[A.obs["pt_bin3"] == "P1"].copy()
    sym2gid = {}
    for gid, s in zip(A.var_names.astype(str), A.var["sym"].astype(str)):
        sym2gid.setdefault(s, gid)
    genes = [("ccnd2", "ccnd2"), ("LOC113823945", "hemocyte TGM")]
    positions = [(0.82, 1.18), (1.82, 2.18)]
    rng = np.random.default_rng(2)
    for (sym, lab), (p0, p1) in zip(genes, positions):
        gid = sym2gid[sym]
        v = early[:, gid].X.toarray().ravel().astype(float)
        vals = [v[(early.obs["group"].astype(str) == g).values] for g in ("blank", "treated")]
        vp = ax.violinplot(vals, positions=[p0, p1], widths=0.30, showextrema=False, showmeans=False)
        for body, col in zip(vp["bodies"], (CTRL, TRT)):
            body.set_facecolor(col); body.set_edgecolor("none"); body.set_alpha(0.55)
        for pos, val, col in zip((p0, p1), vals, (CTRL, TRT)):
            jitter = rng.uniform(-0.045, 0.045, size=len(val))
            ax.scatter(pos + jitter, val, s=1.0, c=col, alpha=0.42, linewidths=0, rasterized=True)
            ax.scatter([pos], [np.median(val)], s=11, c="white", edgecolors=col, linewidths=0.6, zorder=5)
        p = mannwhitneyu(vals[1], vals[0], alternative="two-sided").pvalue
        ymax = max(np.max(vals[0]), np.max(vals[1])) + 0.55
        ax.plot([p0, p1], [ymax, ymax], color="#333333", lw=0.5)
        ax.text((p0 + p1) / 2, ymax + 0.12, "P = %.1e" % p, ha="center", va="bottom", fontsize=5.0)
    ax.set_xticks([1.0, 2.0]); ax.set_xticklabels(["ccnd2", "hemocyte TGM"], fontsize=5.4)
    ax.set_ylabel("Log-normalized expression", fontsize=5.6)
    ax.set_xlim(0.55, 2.45); ax.set_ylim(-0.1, 6.2)
    ax.tick_params(axis="y", labelsize=5.0, length=1.2); ax.tick_params(axis="x", length=0)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.set_title("Early PH: matched pseudotime bin", fontsize=6.6, pad=4)
    ax.legend(handles=[Patch(color=CTRL, label="Control"), Patch(color=TRT, label="Hemospasia")],
              loc="upper left", frameon=False, fontsize=5.0, handlelength=0.8,
              handletextpad=0.3, borderpad=0.1)


def box(ax, x, y, w, h, text, fc, ec, fs=5.2, tc="white"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                       linewidth=0.6, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=3, linespacing=1.05)


def arrow(ax, p1, p2, color="#333333", lw=0.9):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>,head_length=0.5,head_width=0.22",
                                 mutation_scale=9, lw=lw, color=color, zorder=4,
                                 shrinkA=1, shrinkB=1))


def draw_model(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    # upstream module
    ax.add_patch(FancyBboxPatch((0.18, 0.82), 0.44, 0.09, boxstyle="round,pad=0.01,rounding_size=0.02",
                                lw=0.6, ec="#888888", fc="#F4F4F4", ls=(0, (3, 2)), zorder=1))
    ax.text(0.40, 0.865, "Candidate upstream: br / Trl / Cfp1 / Taf7 / Kdm4B",
            ha="center", va="center", fontsize=4.9, color="#333333", zorder=2)
    box(ax, 0.02, 0.31, 0.16, 0.22, "Hemospasia", "#BDBDBD", "#777777", fs=5.0, tc="#222222")
    box(ax, 0.22, 0.31, 0.15, 0.22, "Early PH", "#B2182B", "#8C1A1A", fs=5.2)
    box(ax, 0.44, 0.57, 0.17, 0.18, "TGM ↑\nprogenitor identity", "#D6604D", "#A33A32", fs=4.8)
    box(ax, 0.44, 0.09, 0.17, 0.18, "ccnd2 ↑\nG1/S entry", "#D6604D", "#A33A32", fs=4.8)
    box(ax, 0.70, 0.31, 0.17, 0.22, "PH expansion", "#4C78A8", "#2E5F85", fs=5.2)
    arrow(ax, (0.18, 0.42), (0.22, 0.42))
    arrow(ax, (0.37, 0.45), (0.44, 0.64), color="#A33A32")
    arrow(ax, (0.37, 0.39), (0.44, 0.20), color="#A33A32")
    arrow(ax, (0.61, 0.66), (0.70, 0.47), color="#A33A32")
    arrow(ax, (0.61, 0.18), (0.70, 0.37), color="#A33A32")
    # negative / differentiation arm, relieved in treated PH
    box(ax, 0.70, 0.015, 0.20, 0.13,
        "Differentiation / restraining arm\nTHBS4↓  PEN-3↓  crustin↓  WAP↓",
        "#E8F0F8", "#4575B4", fs=4.7, tc="#1F4E79")
    # blunt inhibitory connection to PH expansion
    ax.plot([0.80, 0.80], [0.145, 0.305], color="#4575B4", lw=0.7, ls=(0, (2, 2)), zorder=4)
    ax.plot([0.785, 0.815], [0.305, 0.305], color="#4575B4", lw=1.2, zorder=5)
    ax.text(0.815, 0.215, "relieved", fontsize=4.6, color="#4575B4", ha="left", va="center")


if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 4.8), facecolor="white")
    axa = fig.add_axes([0.07, 0.40, 0.40, 0.50])
    axb = fig.add_axes([0.55, 0.40, 0.40, 0.50])
    axc = fig.add_axes([0.035, 0.035, 0.94, 0.28])
    draw_volcano(axa)
    draw_violin(axb)
    draw_model(axc)
    axa.text(-0.18, 1.04, "a", transform=axa.transAxes, fontsize=8, fontweight="bold", va="top")
    axb.text(-0.18, 1.04, "b", transform=axb.transAxes, fontsize=8, fontweight="bold", va="top")
    axc.text(-0.02, 1.02, "c", transform=axc.transAxes, fontsize=8, fontweight="bold", va="top")
    fig.text(0.50, 0.012,
             "Candidate model; the differentiation/restraining arm is suppressed, not a proven direct interaction.",
             ha="center", va="bottom", fontsize=4.8, color="#555555")
    p = os.path.join(PREV, "04b_ccnd2_TGM_mechanism_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)
