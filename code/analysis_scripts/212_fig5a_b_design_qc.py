# -*- coding: utf-8 -*-
"""212_fig5a_b_design_qc.py
Draw only the two missing panels of the first single-cell summary figure:
a) study design and cell yield
b) QC distributions of retained cells

Existing t-SNE, heatmap, composition and cell-cycle panels are intentionally not redrawn.
Exports standalone editable SVG + 1200-dpi TIFF panels and one combined review preview.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image

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

CTRL = "#4C78A8"
TRT = "#D6604D"
CTRL_LIGHT = "#EAF2F8"
TRT_LIGHT = "#FBEAE8"
INK = "#222222"
GREY = "#666666"
CELL_COL = ["#B2182B", "#E08214", "#D6604D", "#8C510A", "#7B3294",
            "#006D2C", "#41AB5D", "#78C679", "#92C5DE", "#2166AC"]


def flatten_tiff(path, dpi=1200):
    with Image.open(path) as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA")
            bg = Image.new("RGB", rgba.size, (255, 255, 255))
            bg.paste(rgba, mask=rgba.getchannel("A"))
            tmp = path + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(dpi, dpi), compression="tiff_lzw")
            os.replace(tmp, path)


def save_panel(fig, stem, preview_dpi=400):
    svg = stem + ".svg"
    tif = stem + ".tiff"
    png = stem + "_preview.png"
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    fig.savefig(stem + "_qa.pdf", bbox_inches="tight", pad_inches=0.03,
                facecolor="white")
    fig.savefig(tif, dpi=1200, bbox_inches="tight", pad_inches=0.03,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    flatten_tiff(tif, 1200)
    fig.savefig(png, dpi=preview_dpi, bbox_inches="tight", pad_inches=0.03,
                facecolor="white")
    print("saved", svg)
    print("saved", tif)
    print("saved", png)


def load_counts():
    raw = pd.read_csv(os.path.join(OUT, "qc_per_cell_raw.csv"))
    A = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"), backed="r")
    obs = A.obs.copy()
    A.file.close()
    n_raw = raw["sample"].value_counts().to_dict()
    n_ret = obs["group"].astype(str).value_counts().to_dict()
    rows = [
        {"group": "Control", "time_h": 0, "n_profiled": int(n_raw.get("blank", 0)),
         "n_qc_retained": int(n_ret.get("blank", 0))},
        {"group": "Hemospasia", "time_h": 6, "n_profiled": int(n_raw.get("treated", 0)),
         "n_qc_retained": int(n_ret.get("treated", 0))},
    ]
    D = pd.DataFrame(rows)
    D["retention_pct"] = 100 * D["n_qc_retained"] / D["n_profiled"]
    D.to_csv(os.path.join(OUT, "Fig5a_study_design_counts.csv"), index=False, encoding="utf-8-sig")
    return D, obs


def load_qc(obs):
    # Use the original pre-SoupX QC values for the retained barcodes.
    q = pd.read_csv(os.path.join(OUT, "qc_per_cell_raw.csv")).set_index("Cells")
    q = q.loc[obs.index].copy()
    q["group"] = obs["group"].astype(str).map({"blank": "Control", "treated": "Hemospasia"}).values
    q = q.rename(columns={"nCount": "nCount", "nFeature": "nFeature", "pct_mt": "pct_mt"})
    n_before = len(q)
    q = q.dropna(subset=["nFeature", "nCount", "pct_mt"]).copy()
    n_after = len(q)
    print(f"QC rows: before={n_before}, after={n_after}, excluded={n_before-n_after}")
    rows = []
    for g, sub in q.groupby("group"):
        row = {"group": g, "n": int(len(sub))}
        for m in ("nFeature", "nCount", "pct_mt"):
            row[m + "_median"] = float(sub[m].median())
            row[m + "_q25"] = float(sub[m].quantile(0.25))
            row[m + "_q75"] = float(sub[m].quantile(0.75))
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "Fig5b_qc_summary.csv"), index=False,
                              encoding="utf-8-sig")
    return q


def draw_panel_a(ax, D):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(-0.015, 1.04, "a", transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")

    # Sample source
    ax.text(0.50, 0.955, "Hemolymph scRNA-seq", ha="center", va="top",
            fontsize=6.6, fontweight="bold", color=INK)
    ax.text(0.50, 0.905, "Litopenaeus vannamei", ha="center", va="top",
            fontsize=5.3, color=GREY)

    # Split line
    ax.plot([0.50, 0.50], [0.865, 0.795], color="#777777", lw=0.7)
    ax.plot([0.27, 0.73], [0.795, 0.795], color="#777777", lw=0.7)
    ax.plot([0.27, 0.27], [0.795, 0.745], color="#777777", lw=0.7)
    ax.plot([0.73, 0.73], [0.795, 0.745], color="#777777", lw=0.7)

    boxes = [
        (0.27, "Control", "0 h", D.iloc[0], CTRL, CTRL_LIGHT),
        (0.73, "Hemospasia", "6 h", D.iloc[1], TRT, TRT_LIGHT),
    ]
    for x, name, time, row, col, face in boxes:
        box = FancyBboxPatch((x - 0.185, 0.585), 0.37, 0.155,
                             boxstyle="round,pad=0.008,rounding_size=0.018",
                             facecolor=face, edgecolor=col, linewidth=0.8)
        ax.add_patch(box)
        ax.text(x, 0.705, f"{name} ({time})", ha="center", va="center",
                fontsize=5.9, fontweight="bold", color=col)
        ax.text(x, 0.655, f"profiled  {row.n_profiled:,}", ha="center", va="center",
                fontsize=5.0, color=INK)
        ax.text(x, 0.612, f"QC-retained  {row.n_qc_retained:,}", ha="center", va="center",
                fontsize=5.2, color=INK)
        ax.add_patch(FancyArrowPatch((x, 0.582), (x, 0.380),
                                     arrowstyle="-|>", mutation_scale=6.0,
                                     lw=0.75, color=col, shrinkA=0, shrinkB=0))

    # Central analysis box
    central = FancyBboxPatch((0.135, 0.105), 0.73, 0.275,
                             boxstyle="round,pad=0.010,rounding_size=0.020",
                             facecolor="#F7F7F7", edgecolor="#777777", linewidth=0.8)
    ax.add_patch(central)
    ax.text(0.50, 0.345, "10× scRNA-seq  ·  9,658 QC-retained cells", ha="center",
            va="top", fontsize=6.0, fontweight="bold", color=INK)
    xs = np.linspace(0.35, 0.65, 10)
    for x, c in zip(xs, CELL_COL):
        ax.scatter([x], [0.290], s=10, c=c, edgecolors="white", linewidths=0.25, zorder=3)
    ax.text(0.50, 0.245, "10 transcriptionally distinct cell states", ha="center",
            va="center", fontsize=5.4, color=INK)
    ax.text(0.50, 0.190, "markers  ·  t-SNE  ·  composition  ·  cell cycle", ha="center",
            va="center", fontsize=5.0, color=GREY)
    ax.text(0.50, 0.145, "QC thresholds selected per cell state", ha="center",
            va="center", fontsize=5.0, color=GREY)

    total_raw = int(D.n_profiled.sum()); total_ret = int(D.n_qc_retained.sum())
    ax.text(0.50, 0.055, f"{total_raw:,} profiled → {total_ret:,} retained "
                          f"({100 * total_ret / total_raw:.1f}%)",
            ha="center", va="center", fontsize=5.2, color=INK)


def draw_panel_b(axes, q):
    metrics = [("nFeature", "Genes/cell"), ("nCount", "UMIs/cell"), ("pct_mt", "Mito. UMI (%)")]
    groups = ["Control", "Hemospasia"]
    colors = [CTRL, TRT]
    axes[0].text(-0.16, 1.14, "b", transform=axes[0].transAxes, fontsize=9,
                 fontweight="bold", va="top", ha="left")
    axes[0].set_title("QC metrics of retained cells  (n = 9,658)", fontsize=6.4,
                      pad=5, weight="bold")
    for i, (ax, (metric, label)) in enumerate(zip(axes, metrics)):
        vals = [q.loc[q.group == g, metric].dropna().to_numpy(float) for g in groups]
        vp = ax.violinplot(vals, positions=[0, 1], widths=0.62, showextrema=False,
                           showmeans=False)
        for body, col in zip(vp["bodies"], colors):
            body.set_facecolor(col)
            body.set_edgecolor(col)
            body.set_alpha(0.55)
            body.set_linewidth(0.4)
        for pos, v, col in zip([0, 1], vals, colors):
            med = float(np.median(v))
            ax.hlines(med, pos - 0.18, pos + 0.18, color="white", lw=2.0, zorder=4)
            ax.hlines(med, pos - 0.18, pos + 0.18, color=INK, lw=0.65, zorder=5)
        lo = float(np.nanpercentile(np.concatenate(vals), 1))
        hi = float(np.nanpercentile(np.concatenate(vals), 99))
        pad = 0.06 * (hi - lo if hi > lo else 1)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(-0.58, 1.58)
        ax.set_ylabel(label, fontsize=5.3, labelpad=2)
        ax.tick_params(axis="y", labelsize=5.0, length=1.5, width=0.45)
        ax.tick_params(axis="x", length=0)
        if i < len(axes) - 1:
            ax.set_xticks([0, 1]); ax.set_xticklabels([])
        else:
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Control\nn = 5,349", "Hemospasia\nn = 4,309"],
                               fontsize=5.0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_linewidth(0.55)
        ax.spines["bottom"].set_linewidth(0.55)
        ax.grid(axis="y", color="#E8E8E8", lw=0.35, zorder=0)
        ax.set_axisbelow(True)
    axes[-1].text(0.50, -0.50,
                  "Main QC: pct_mt <10%, nFeature ≥200; no doublet filtering.\n"
                  "Values before SoupX; violins 1st–99th percentile; line = median.",
                  transform=axes[-1].transAxes, ha="center", va="top", fontsize=5.0,
                  color=GREY)


def make_panel_a(D, path):
    fig = plt.figure(figsize=(3.55, 2.85), facecolor="white")
    ax = fig.add_axes([0.025, 0.035, 0.95, 0.90])
    draw_panel_a(ax, D)
    save_panel(fig, path)
    plt.close(fig)


def make_panel_b(q, path):
    fig = plt.figure(figsize=(3.75, 2.85), facecolor="white")
    gs = fig.add_gridspec(3, 1, left=0.16, right=0.975, top=0.865, bottom=0.165,
                         hspace=0.12)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    draw_panel_b(axes, q)
    save_panel(fig, path)
    plt.close(fig)


def make_combined(D, q, path):
    fig = plt.figure(figsize=(7.20, 2.95), facecolor="white")
    gs = fig.add_gridspec(1, 2, left=0.035, right=0.985, top=0.955, bottom=0.075,
                          width_ratios=[1.0, 1.18], wspace=0.12)
    axa = fig.add_subplot(gs[0, 0])
    gsb = gs[0, 1].subgridspec(3, 1, hspace=0.12)
    axes_b = [fig.add_subplot(gsb[i, 0]) for i in range(3)]
    draw_panel_a(axa, D)
    draw_panel_b(axes_b, q)
    fig.savefig(path, dpi=450, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    print("saved", path)
    plt.close(fig)


def main():
    D, obs = load_counts()
    q = load_qc(obs)
    print(D.to_string(index=False))
    print(q.groupby("group")[["nFeature", "nCount", "pct_mt"]].median().to_string())
    make_panel_a(D, os.path.join(PARTS, "Fig5a_study_design"))
    make_panel_b(q, os.path.join(PARTS, "Fig5b_qc"))
    make_combined(D, q, os.path.join(PARTS, "Fig5ab_missing_panels_preview.png"))


if __name__ == "__main__":
    main()
