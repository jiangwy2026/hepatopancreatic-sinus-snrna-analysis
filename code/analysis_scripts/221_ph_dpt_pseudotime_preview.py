# -*- coding: utf-8 -*-
"""221_ph_dpt_pseudotime_preview.py
PH-specific DPT pseudotime previews:
1) density of DPT pseudotime by PH substate and group;
2) conventional diffusion-map pseudotime plot with PH substate labels and root.
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
from scipy.stats import gaussian_kde, mannwhitneyu
from matplotlib.lines import Line2D
from PIL import Image

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-20", "PH亚亚群拟时")
os.makedirs(PARTS, exist_ok=True); os.makedirs(DELIV, exist_ok=True)
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

PH_KEYS = ["PH_s0", "PH_s1", "PH_s2", "PH_s3"]
PH_LABELS = {"PH_s0": "s0", "PH_s1": "s1",
             "PH_s2": "s2", "PH_s3": "s3"}
COL = {"PH_s0": "#59A14F", "PH_s1": "#E15759", "PH_s2": "#F28E2B", "PH_s3": "#4E79A7"}


def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan)
    ok = np.isfinite(p); idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]; m = len(order); prev = 1.0
    for rank, i in enumerate(order[::-1], 1):
        val = min(prev, p[i] * m / (m - rank + 1)); out[i] = val; prev = val
    return out


def flatten_tiff(path, dpi=1200):
    with Image.open(path) as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA")
            bg = Image.new("RGB", rgba.size, (255, 255, 255))
            bg.paste(rgba, mask=rgba.getchannel("A"))
            tmp = path + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(dpi, dpi), compression="tiff_lzw")
            os.replace(tmp, path)


def save_export(fig, name):
    png = os.path.join(PARTS, name + "_preview.png")
    svg_local = os.path.join(PARTS, name + ".svg")
    tif_local = os.path.join(PARTS, name + ".tiff")
    svg_deliv = os.path.join(DELIV, name + ".svg")
    tif_deliv = os.path.join(DELIV, name + ".tiff")
    fig.savefig(png, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(svg_local, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(tif_local, dpi=1200, bbox_inches="tight", pad_inches=0.04,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    flatten_tiff(tif_local, 1200)
    fig.savefig(svg_deliv, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(tif_deliv, dpi=1200, bbox_inches="tight", pad_inches=0.04,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    flatten_tiff(tif_deliv, 1200)
    print("saved", png)
    print("saved", svg_local)
    print("saved", tif_local)
    print("saved", svg_deliv)
    print("saved", tif_deliv)


def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    sub = pd.read_csv(os.path.join(OUT, "dynamic_subpop_subclusters", "PH_subcluster_umap.csv"), index_col=0)
    sub["sub"] = "PH_s" + sub["subcluster"].astype(str).str.replace(".0", "", regex=False)
    obs = A.obs.join(sub[["sub"]], how="left")
    obs["group2"] = obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    ph = obs[obs["ct"].astype(str) == "PH"].copy()
    pt = ph["dpt_pseudotime"].astype(float)
    print("PH DPT cells", len(ph), "range", float(pt.min()), float(pt.max()))

    pvals = []
    for key in PH_KEYS:
        c = pt[(ph["sub"].astype(str) == key) & (ph.group2 == "Control")].dropna()
        t = pt[(ph["sub"].astype(str) == key) & (ph.group2 == "Hemospasia")].dropna()
        pvals.append(mannwhitneyu(c, t, alternative="two-sided").pvalue)
    fdr = bh(pvals)
    sig = {PH_KEYS[i] for i in range(len(PH_KEYS)) if fdr[i] < 0.05}
    summ = []
    for key in PH_KEYS:
        for g in ("Control", "Hemospasia"):
            x = pt[(ph["sub"].astype(str) == key) & (ph.group2 == g)].dropna()
            summ.append({"key": key, "group": g, "n": len(x), "median_dpt": float(x.median()),
                         "q25": float(x.quantile(.25)), "q75": float(x.quantile(.75)),
                         "mannwhitney_p": pvals[PH_KEYS.index(key)], "fdr": fdr[PH_KEYS.index(key)]})
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "Fig5_PH_dpt_pseudotime_summary.csv"),
                              index=False, encoding="utf-8-sig")

    # Density plot
    fig = plt.figure(figsize=(4.9, 3.8), facecolor="white")
    ax = fig.add_axes([0.27, 0.235, 0.68, 0.635])
    xmin, xmax = pt.quantile([0.01, 0.99]).astype(float).values
    pad = 0.12 * (xmax - xmin); xx = np.linspace(xmin - pad, xmax + pad, 240)
    for row, key in enumerate(PH_KEYS[::-1]):
        d = ph[ph["sub"].astype(str) == key]
        for g, ls, alpha in [("Control", (0, (3, 2)), 0.80), ("Hemospasia", "-", 1.0)]:
            x = d.loc[d.group2 == g, "dpt_pseudotime"].dropna().astype(float).values
            if len(x) >= 5:
                y = gaussian_kde(x, bw_method="scott")(xx); y = y / max(y.max(), 1e-9) * 0.80
                ax.plot(xx, row + y, color=COL[key], lw=0.95, ls=ls, alpha=alpha)
                ax.fill_between(xx, row, row + y, color=COL[key], alpha=0.08 if g == "Hemospasia" else 0)
                ax.plot([np.median(x), np.median(x)], [row - 0.05, row + 0.10], color=COL[key],
                        lw=1.0 if g == "Hemospasia" else 0.7, ls=ls)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([PH_LABELS[k] + (" *" if k in sig else "") for k in PH_KEYS[::-1]], fontsize=5.2)
    ax.set_xlim(xx.min(), xx.max()); ax.set_ylim(-0.35, 3.95)
    ax.set_xlabel("DPT pseudotime", fontsize=5.8)
    ax.tick_params(axis="x", labelsize=5.1, length=1.4, width=0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#EFEFEF", lw=0.35, zorder=0)
    ax.set_title("PH substate DPT pseudotime", fontsize=7.0, pad=5)
    handles = [Line2D([0], [0], color="#444444", lw=0.9, label="Hemospasia"),
               Line2D([0], [0], color="#444444", lw=0.9, ls=(0, (3, 2)), label="Control"),
               Line2D([0], [0], color="#444444", lw=0, marker="|", markersize=8, label="Median")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.085), ncol=3,
               frameon=False, fontsize=5.0, handlelength=1.2, columnspacing=1.0,
               handletextpad=0.35)
    fig.text(0.5, 0.018, "* FDR < 0.05 for group difference (cell-level Mann–Whitney). "
             "DPT pseudotime is a transcriptional ordering, not lineage tracing.",
             ha="center", va="bottom", fontsize=4.9, color="#555555")
    save_export(fig, "Fig5_PH_DPT_pseudotime_density")
    plt.close(fig)

    # Conventional diffusion-map pseudotime plot
    U = np.asarray(A.obsm["X_diffmap"])[:, :2]
    keep = np.where(A.obs.index.isin(ph.index))[0]
    U = U[keep]
    order = pd.Index(A.obs_names[keep])
    d = ph.loc[order].copy()
    root_i = int(np.nanargmin(d["dpt_pseudotime"].values))
    fig = plt.figure(figsize=(4.6, 4.0), facecolor="white")
    ax = fig.add_axes([0.07, 0.10, 0.76, 0.80])
    pt_rank = d["dpt_pseudotime"].astype(float).rank(pct=True).values
    sca = ax.scatter(U[:, 0], U[:, 1], c=pt_rank, s=3.2,
                     cmap="viridis", linewidths=0, rasterized=True, alpha=0.90,
                     vmin=0, vmax=1)
    for key in PH_KEYS:
        m = d["sub"].astype(str).values == key
        if m.sum() >= 3:
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), key.split("_")[1],
                    fontsize=5.2, ha="center", va="center", color="white", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.12", fc="#333333", ec="none", alpha=0.72), zorder=8)
    ax.scatter([U[root_i, 0]], [U[root_i, 1]], marker="*", s=95, c="#FFD92F",
               edgecolors="#222222", linewidths=0.5, zorder=10)
    ax.annotate("root", xy=(U[root_i, 0], U[root_i, 1]), xytext=(5, -9),
                textcoords="offset points", fontsize=5.0, color="#222222", ha="left", va="top")
    x1, x2 = np.quantile(U[:, 0], [0.01, 0.99]); y1, y2 = np.quantile(U[:, 1], [0.01, 0.99])
    xp = 0.08 * (x2 - x1); yp = 0.08 * (y2 - y1)
    ax.set_xlim(x1 - xp, x2 + xp); ax.set_ylim(y1 - yp, y2 + yp)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.4); sp.set_color("#777777")
    ax.set_title("PH diffusion pseudotime", fontsize=7.0, pad=5)
    cax = fig.add_axes([0.86, 0.20, 0.020, 0.58])
    cb = fig.colorbar(sca, cax=cax); cb.set_label("Within-PH DPT percentile", fontsize=5.5, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.3, width=0.4); cb.outline.set_linewidth(0.4)
    fig.text(0.5, 0.02, "Diffusion map from the existing hemocyte trajectory object; colour = within-PH DPT percentile. "
             "Root = earliest PH cell.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out2 = os.path.join(PARTS, "Fig5_PH_dpt_pseudotime_diffmap_preview.png")
    fig.savefig(out2, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out2)
    plt.close(fig)


if __name__ == "__main__":
    main()
