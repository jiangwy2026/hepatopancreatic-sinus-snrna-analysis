# -*- coding: utf-8 -*-
"""176_sc_summary_new_names_preview.py
Four-panel single-cell summary using the final subtype names:
a composition, b PH-rooted pseudotime, c velocity summary, d ccnd2/TGM mechanism.
Preview only; final SVG/TIFF after user confirmation.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import rankdata

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names")
os.makedirs(PREV, exist_ok=True)

COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}
ORDER = ["PH", "GH", "TS", "MH2"]
CTRL = "#4C78A8"; TRT = "#D6604D"; GREY = "#8A8A8A"
MAP = {"MH1": "TS", "PH": "PH", "GH": "GH", "MH2": "MH2"}


def clean_axes(ax, xlabel=None, ylabel=None):
    ax.set_xticks([]); ax.set_yticks([])
    if xlabel: ax.set_xlabel(xlabel, fontsize=6.0)
    if ylabel: ax.set_ylabel(ylabel, fontsize=6.0)
    for sp in ax.spines.values():
        sp.set_color("#777777"); sp.set_linewidth(0.45)


def normalize_rank(x):
    x = np.asarray(x, dtype=float)
    good = np.isfinite(x)
    out = np.full(len(x), np.nan)
    if good.sum():
        r = rankdata(x[good], method="average")
        out[good] = (r - 1) / max(1, good.sum() - 1)
    return out


def draw_composition(ax):
    comp = pd.read_csv(os.path.join(OUT, "composition_new_names_10groups.csv"))
    chi = json.load(open(os.path.join(OUT, "composition_new_names_chi2.json")))
    y = np.arange(len(comp))
    ax.barh(y - 0.16, comp.pct_blank, height=0.32, color=CTRL, edgecolor="none", label="Control (0 h)")
    ax.barh(y + 0.16, comp.pct_treated, height=0.32, color=TRT, edgecolor="none", label="Hemospasia (6 h)")
    labels = []
    for _, r in comp.iterrows():
        lab = "%s (%s)" % (r.group, format(int(r.n_total), ","))
        if r.group in ("PH", "MH2", "PM-like", "R"):
            lab += "; %.2f×" % r.fold_treated_vs_blank
        labels.append(lab)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.7)
    ax.invert_yaxis()
    ax.set_xlim(0, 55)
    ax.set_xlabel("Cells (%%)   (χ² = %.1f; $P$ = %.1f × 10$^{-135}$)" % (chi["chi2"], chi["p"] / 1e-135),
                  fontsize=5.5)
    ax.set_title("Cell composition (%)", fontsize=7.0, pad=5)
    ax.legend(loc="lower right", frameon=False, fontsize=5.3, handlelength=0.9,
              handletextpad=0.35, borderpad=0.1, labelspacing=0.35)
    ax.tick_params(axis="x", labelsize=5.3, length=1.7, width=0.5)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#E5E5E5", lw=0.35, zorder=0)
    ax.set_axisbelow(True)


def draw_pseudotime(ax):
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    A.obs["ct_new"] = pd.Categorical(A.obs["ct"].astype(str).map(MAP), categories=ORDER)
    ct = A.obs["ct_new"].astype(str).values
    U = np.asarray(A.obsm["X_umap"])
    root = int(A.uns.get("iroot", 0))
    if ct[root] != "PH":
        ph = np.where(ct == "PH")[0]
        root = int(ph[np.argmin(((U[ph] - U[ph].mean(0)) ** 2).sum(1))])
    dpt = np.asarray(A.obs["dpt_pseudotime"], dtype=float)
    dpt[~np.isfinite(dpt)] = np.nan
    pt = normalize_rank(dpt)
    o = np.argsort(-pt)
    s = ax.scatter(U[o, 0], U[o, 1], c=pt[o], s=1.15, cmap="viridis", linewidths=0,
                   rasterized=True, vmin=0, vmax=1)
    cent = {g: U[ct == g].mean(0) for g in ORDER}
    # Candidate paths based on DPT ordering. These are not velocity directions.
    paths = [(cent["PH"], cent["TS"], -0.10, "#2B2B2B", 1.0),
             (cent["TS"], cent["GH"], 0.30, "#2166AC", 0.85),
             (cent["TS"], cent["MH2"], -0.30, "#8C510A", 0.85)]
    for p1, p2, rad, col, lw in paths:
        ax.add_patch(FancyArrowPatch(p1, p2, connectionstyle="arc3,rad=%s" % rad,
                                     arrowstyle="-|>,head_length=0.45,head_width=0.20",
                                     mutation_scale=8, lw=lw, color=col, alpha=0.9,
                                     zorder=9, shrinkA=2, shrinkB=2,
                                     path_effects=[pe.withStroke(linewidth=2.0, foreground="white")]))
    ax.scatter(U[root, 0], U[root, 1], s=58, marker="*", c="#FF3B30", edgecolors="white",
               linewidths=0.6, zorder=11)
    ax.annotate("root (PH)", xy=(U[root, 0], U[root, 1]), xytext=(0, 13), textcoords="offset points",
                fontsize=5.8, color="#B2182B", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.55, color="#B2182B"), zorder=12)
    offs = {"PH": (-12, 9), "GH": (0, 0), "TS": (0, 0), "MH2": (0, 0)}
    for g in ORDER:
        c = cent[g]
        ax.annotate(g, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=6.0, weight="bold", ha="center", va="center", color="white",
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="#333333")], zorder=10)
    cb = fig.colorbar(s, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("DPT pseudotime (rank display)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5.0, length=1.0); cb.outline.set_visible(False)
    clean_axes(ax, "UMAP 1", "UMAP 2")
    ax.set_title("PH-rooted pseudotime", fontsize=7.0, pad=5)
    ax.text(-0.17, 1.045, "b", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")


def draw_velocity(ax_heat, ax_pt):
    V = pd.read_csv(os.path.join(OUT, "velocity_flow_4groups_offdiag.csv"), index_col=0)
    V = V.loc[ORDER, ORDER]
    A = V.to_numpy(copy=True)
    M = np.ma.masked_invalid(A)
    cmap = plt.cm.magma.copy(); cmap.set_bad("#F2F2F2")
    im = ax_heat.imshow(M, cmap=cmap, vmin=0, vmax=0.37, aspect="auto", interpolation="nearest")
    ax_heat.set_xticks(range(4)); ax_heat.set_yticks(range(4))
    ax_heat.set_xticklabels(ORDER, fontsize=5.5); ax_heat.set_yticklabels(ORDER, fontsize=5.5)
    ax_heat.xaxis.tick_top(); ax_heat.tick_params(length=0)
    ax_heat.set_xlabel("Target", fontsize=5.6, labelpad=4); ax_heat.set_ylabel("Source", fontsize=5.6)
    for i in range(4):
        for j in range(4):
            if i == j or not np.isfinite(A[i, j]):
                continue
            ax_heat.text(j, i, "%.1f" % (100 * A[i, j]), ha="center", va="center",
                         fontsize=5.2, color="white" if A[i, j] < 0.25 else "black")
    cb = fig.colorbar(im, ax=ax_heat, fraction=0.045, pad=0.03)
    cb.set_label("Fraction of source transitions", fontsize=5.3)
    cb.ax.tick_params(labelsize=5.0, length=1.0); cb.outline.set_visible(False)
    ax_heat.set_title("Velocity flow (self excluded)", fontsize=7.0, pad=18)
    ax_heat.text(0.0, -0.24, "PH→TS 30.1% vs TS→PH 5.3% (5.7×)",
                 transform=ax_heat.transAxes, fontsize=5.3, color="#B2182B")

    VP = pd.read_csv(os.path.join(OUT, "velocity_per_group_new_names.csv"))
    VP = VP.set_index("group_new").loc[ORDER]
    yy = np.arange(len(ORDER))
    colors = [COL[g] for g in ORDER]
    ax_pt.barh(yy, VP["velocity_pseudotime"], color=colors, height=0.55, edgecolor="none")
    ax_pt.set_yticks(yy); ax_pt.set_yticklabels(ORDER, fontsize=5.4)
    ax_pt.invert_yaxis()
    ax_pt.set_xlim(0.105, 0.145)
    ax_pt.set_xlabel("Velocity pseudotime", fontsize=5.7)
    ax_pt.tick_params(axis="x", labelsize=5.0, length=1.2)
    ax_pt.tick_params(axis="y", length=0)
    for i, v in enumerate(VP["velocity_pseudotime"]):
        ax_pt.text(v + 0.001, i, "%.3f" % v, va="center", fontsize=5.0)
    ax_pt.set_title("Velocity pseudotime", fontsize=5.7, pad=3)
    for sp in ("top", "right"):
        ax_pt.spines[sp].set_visible(False)
    ax_pt.text(-0.27, 1.45, "c", transform=ax_heat.transAxes, fontsize=8, fontweight="bold", va="top")


def draw_mechanism(ax_vol, ax_viol):
    D = pd.read_csv(os.path.join(OUT, "PH_early_treated_vs_blank_DE.csv"))
    D = D[np.isfinite(D.logfoldchanges) & np.isfinite(D.pvals_adj)].copy()
    D["neglogq"] = -np.log10(np.clip(D.pvals_adj.values, 1e-300, 1))
    sig = (D.pvals_adj < 0.05) & (D.logfoldchanges.abs() > 0.5)
    up = sig & (D.logfoldchanges > 0); dn = sig & (D.logfoldchanges < 0)
    ax_vol.scatter(D.loc[~sig, "logfoldchanges"], D.loc[~sig, "neglogq"], s=2.5, c="#BDBDBD",
                   linewidths=0, alpha=0.55, rasterized=True)
    ax_vol.scatter(D.loc[up, "logfoldchanges"], D.loc[up, "neglogq"], s=5.0, c="#D73027",
                   linewidths=0, alpha=0.85, rasterized=True)
    ax_vol.scatter(D.loc[dn, "logfoldchanges"], D.loc[dn, "neglogq"], s=5.0, c="#4575B4",
                   linewidths=0, alpha=0.75, rasterized=True)
    ax_vol.axvline(-0.5, color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax_vol.axvline(0.5, color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax_vol.axhline(-np.log10(0.05), color="#777777", lw=0.4, ls=(0, (2, 2)))
    ax_vol.set_xlim(-6.5, 3.5); ax_vol.set_ylim(0, 30)
    ax_vol.set_xlabel("log$_2$ fold change (treated / control)", fontsize=5.7)
    ax_vol.set_ylabel("−log$_{10}$ adjusted $P$", fontsize=5.7)
    ax_vol.tick_params(labelsize=5.0, length=1.2)
    for sp in ("top", "right"): ax_vol.spines[sp].set_visible(False)
    ax_vol.set_title("Early PH response: treated vs control", fontsize=7.0, pad=4)
    labels = {
        "LOC113823945": ("hemocyte TGM", (0.25, 0.5), "#B2182B"),
        "ccnd2": ("ccnd2", (0.25, -0.7), "#B2182B"),
        "SSPO": ("SSPO", (-0.2, 1.2), "#B2182B"),
        "LOC113808957": ("Notch-like", (0.25, 1.8), "#B2182B"),
        "Igfbp7": ("Igfbp7", (-0.3, 1.0), "#B2182B"),
        "THBS4": ("THBS4", (0.2, 1.0), "#2166AC"),
    }
    for gene, (lab, off, col) in labels.items():
        r = D[D.gene == gene]
        if r.empty: continue
        r = r.iloc[0]
        ax_vol.annotate(lab, xy=(r.logfoldchanges, min(r.neglogq, 29.5)),
                        xytext=(r.logfoldchanges + off[0], min(r.neglogq + off[1], 29.5)),
                        fontsize=5.1, color=col, ha="left", va="bottom",
                        arrowprops=dict(arrowstyle="-", lw=0.35, color=col, alpha=0.8))

    # violins for ccnd2 and hemocyte TGM in the early PH bin
    A = ad.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    A.obs["pt_bin3"] = pd.qcut(A.obs["dpt_pseudotime"].astype(float), 3, labels=["P1", "P2", "P3"])
    early = A[A.obs["pt_bin3"] == "P1"].copy()
    sym2gid = {}
    for gid, s in zip(A.var_names.astype(str), A.var["sym"].astype(str)):
        sym2gid.setdefault(s, gid)
    genes = [("ccnd2", "ccnd2"), ("LOC113823945", "hemocyte TGM")]
    positions = [(0.82, 1.18), (1.82, 2.18)]
    rng = np.random.default_rng(1)
    for (sym, lab), (p0, p1) in zip(genes, positions):
        gid = sym2gid[sym]
        v = early[:, gid].X.toarray().ravel().astype(float)
        vals = [v[(early.obs["group"].astype(str) == g).values] for g in ("blank", "treated")]
        vp = ax_viol.violinplot(vals, positions=[p0, p1], widths=0.30, showextrema=False, showmeans=False)
        for body, col in zip(vp["bodies"], (CTRL, TRT)):
            body.set_facecolor(col); body.set_edgecolor("none"); body.set_alpha(0.55)
        for pos, val, col in zip((p0, p1), vals, (CTRL, TRT)):
            jitter = rng.uniform(-0.045, 0.045, size=len(val))
            ax_viol.scatter(pos + jitter, val, s=1.1, c=col, alpha=0.45, linewidths=0, rasterized=True)
            ax_viol.scatter([pos], [np.median(val)], s=10, c="white", edgecolors=col, linewidths=0.6, zorder=5)
        pv = pd.read_csv(os.path.join(OUT, "PH_early_ccnd2_TGM_stats.csv"))
        p = pv[(pv.gene == sym) & (pv.group == "p_value")]["mannwhitney_p"].iloc[0]
        ymax = max(np.max(vals[0]), np.max(vals[1])) + 0.5
        ax_viol.plot([p0, p1], [ymax, ymax], color="#333333", lw=0.5)
        ax_viol.text((p0 + p1) / 2, ymax + 0.12, "P = %.1e" % p, ha="center", va="bottom", fontsize=5.0)
    ax_viol.set_xticks([1.0, 2.0]); ax_viol.set_xticklabels(["ccnd2", "hemocyte TGM"], fontsize=5.4)
    ax_viol.set_ylabel("Log-normalized expression", fontsize=5.6)
    ax_viol.set_xlim(0.55, 2.45)
    ax_viol.tick_params(axis="y", labelsize=5.0, length=1.2)
    ax_viol.tick_params(axis="x", length=0)
    for sp in ("top", "right"): ax_viol.spines[sp].set_visible(False)
    ax_viol.set_title("Early PH: matched pseudotime bin", fontsize=6.0, pad=3)
    ax_viol.legend(handles=[Patch(color=CTRL, label="Control"), Patch(color=TRT, label="Hemospasia")],
                   loc="upper left", frameon=False, fontsize=5.0, handlelength=0.8,
                   handletextpad=0.3, borderpad=0.1)


if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 7.2), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.075, right=0.975, top=0.925, bottom=0.095,
                          wspace=0.38, hspace=0.42, width_ratios=[1.0, 1.05])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axb.set_box_aspect(1)
    gsc = gs[1, 0].subgridspec(2, 1, height_ratios=[2.0, 1.0], hspace=0.55)
    axc1 = fig.add_subplot(gsc[0]); axc2 = fig.add_subplot(gsc[1])
    gsd = gs[1, 1].subgridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.55)
    axd1 = fig.add_subplot(gsd[0]); axd2 = fig.add_subplot(gsd[1])

    draw_composition(axa)
    draw_pseudotime(axb)
    draw_velocity(axc1, axc2)
    draw_mechanism(axd1, axd2)

    axa.text(-0.17, 1.04, "a", transform=axa.transAxes, fontsize=8, fontweight="bold", va="top")
    axd1.text(-0.08, 1.04, "d", transform=axd1.transAxes, fontsize=8, fontweight="bold", va="top")
    fig.text(0.5, 0.018, "Pseudotime arrows indicate DPT-based candidate paths; RNA velocity robustly supports PH as the source.",
             ha="center", va="bottom", fontsize=5.0, color="#333333")
    p = os.path.join(PREV, "Fig_SC_4panel_new_names_preview.png")
    fig.savefig(p, dpi=350, bbox_inches=None, facecolor="white")
    plt.close(fig)
    print("saved", p)

