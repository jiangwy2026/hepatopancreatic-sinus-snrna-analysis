# -*- coding: utf-8 -*-
"""
103_fig2_heatmap_logcells.py - single-cell marker heatmap where the number of cells shown per
cluster is proportional to log10(cluster size), so small clusters stay visible.
Displayed width of a cluster therefore tracks log(cell number), not the raw number.
"""
import sys, io, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
pd.set_option("display.width", 250); sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
C_LOG = 25.0          # cells shown per cluster = round(C_LOG * log10(n)); tune to taste
MIN_SHOW = 20

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B",
            "8": "NE", "9": "Neu"}
    ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
    A.obs["ct"] = pd.Categorical(A.obs.vendor_cluster.astype(str).map(base).replace({"NE": "NEU", "Neu": "NEU"}),
                                 categories=ORD)
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
           "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}
    P = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected.csv"))
    symv = A.var["sym"].astype(str).values; gidv = A.var_names.astype(str)
    sym2i = {}; sym2g = {}
    for i, (s, g) in enumerate(zip(symv, gidv)):
        sym2i.setdefault(s, i); sym2g.setdefault(s, g)
    P["gid"] = P.gene.map(sym2g)

    rng = np.random.default_rng(0)
    cell_idx, marks, shown, total = [], [], [], []
    rows = []
    for g in ORD:
        idx = np.where((A.obs["ct"] == g).values)[0]
        n = len(idx)
        k = max(MIN_SHOW, int(round(C_LOG * math.log10(n))))
        k = min(k, n)
        pick = np.sort(rng.choice(idx, k, replace=False))
        cell_idx += list(pick); marks.append(len(cell_idx)); shown.append(k); total.append(n)
        rows.append({"cluster": g, "n_cells": n, "n_shown": k, "frac": round(k / n, 3)})
    tshow = pd.DataFrame(rows)
    print("--- displayed cells per cluster (log-scaled) ---")
    print(tshow.to_string(index=False))
    print("total cells drawn:", sum(shown), "of", sum(total))
    tshow.to_csv(os.path.join(OUT, "fig2_heatmap_display_counts.csv"), index=False)
    cell_idx = np.array(cell_idx)
    E = A[cell_idx][:, [sym2i[g] for g in P.gene]].X.toarray().astype(np.float32)
    den = np.percentile(E, 99, axis=0) + 1e-9
    Z = np.clip(E / den, 0, 1).T
    print("heatmap matrix: %d genes x %d cells" % Z.shape, flush=True)

    fig = plt.figure(figsize=(9.2, 5.0))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.030, 0.022, 1.0], width_ratios=[0.24, 1.0, 0.026],
                          hspace=0.02, wspace=0.03, left=0.02, right=0.965, top=0.965, bottom=0.03)
    axb = fig.add_subplot(gs[0, 1]); axn = fig.add_subplot(gs[1, 1])
    axh = fig.add_subplot(gs[2, 1]); axl = fig.add_subplot(gs[2, 0]); axc = fig.add_subplot(gs[2, 2])
    nG, nC = Z.shape
    bars = np.zeros((1, nC))
    for i in range(len(ORD)):
        a = 0 if i == 0 else marks[i - 1]; bars[0, a:marks[i]] = i
    axb.imshow(bars, aspect="auto", cmap=ListedColormap([COL[g] for g in ORD]),
               norm=BoundaryNorm(np.arange(-0.5, len(ORD) + 0.5), len(ORD)), extent=(0, nC, 0, 1))
    axb.set_xticks([]); axb.set_yticks([])
    for i, g in enumerate(ORD):
        a = 0 if i == 0 else marks[i - 1]; b = marks[i]
        axb.text((a + b) / 2, 0.5, g, ha="center", va="center", fontsize=5.6, weight="bold",
                 color="white", clip_on=False)
        if i > 0: axb.axvline(a, color="white", lw=1.0)
    for sp in axb.spines.values(): sp.set_visible(False)
    # real cell numbers above the bar
    axn.set_xlim(0, nC); axn.set_ylim(0, 1); axn.axis("off")
    for i, g in enumerate(ORD):
        a = 0 if i == 0 else marks[i - 1]; b = marks[i]
        axn.text((a + b) / 2, 0.0, "n=%s" % format(total[i], ","), ha="center", va="bottom", fontsize=4.2)
    im = axh.imshow(Z, aspect="auto", cmap="magma", vmin=0, vmax=1, interpolation="nearest")
    axh.set_xticks([]); axh.set_yticks([])
    for sp in axh.spines.values(): sp.set_visible(False)
    for i in range(1, len(ORD)):
        axh.axvline(marks[i - 1], color="white", lw=1.0)
    axl.set_xlim(0, 1); axl.set_ylim(nG - 0.5, -0.5); axl.axis("off")
    for i, r in enumerate(P.itertuples()):
        axl.text(1.0, i, "%s (%s)" % (r.gene, r.gid) if pd.notna(r.gid) else r.gene,
                 ha="right", va="center", fontsize=4.8, color="black")
    pos = axh.get_position(); cw, ch = 0.013, pos.height * 0.20
    axc.set_position([pos.x1 + 0.006, pos.y1 - ch - 0.015, cw, ch])
    cb = fig.colorbar(im, cax=axc); cb.set_label("Relative expression", fontsize=5)
    cb.ax.tick_params(labelsize=4.5, length=1.2); cb.outline.set_visible(False)
    bp = os.path.join(FIG, "Fig2_marker_heatmap_singlecell_logcells")
    fig.savefig(bp + "_preview.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    print("saved", bp + "_preview.png")

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\fig2_logcells.log", "w", encoding="utf8", buffering=1)
    main()
