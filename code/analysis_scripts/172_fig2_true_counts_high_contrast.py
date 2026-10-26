# -*- coding: utf-8 -*-
"""172_fig2_true_counts_high_contrast.py
Figure 2 revision:
- clusters ordered from largest to smallest true cell number;
- all 9,658 single cells shown in the main heatmap, with widths proportional to true n;
- MyoEpi (65) and M (45) enlarged at right;
- viridis for stronger contrast.
Exports: editable SVG + 1200-dpi RGB LZW TIFF.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style

apply_pub_style(7)
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
sc.settings.verbosity = 0

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
FIG = os.path.join(OUT, "figures_pub")
PREV = os.path.join(FIG, "final_previews")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-17")
os.makedirs(PREV, exist_ok=True)
os.makedirs(DELIV, exist_ok=True)

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "M": "#2166AC",
       "MyoEpi": "#92C5DE", "PM-like": "#7B3294"}
DPI = 1200


def save_svg_tiff(fig, stem, preview_dpi=320):
    svg = stem + ".svg"
    tif = stem + ".tiff"
    fig.savefig(svg, bbox_inches=None, facecolor="white")
    fig.savefig(tif, dpi=DPI, bbox_inches=None, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tif) as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA")
            bg = Image.new("RGB", rgba.size, (255, 255, 255))
            bg.paste(rgba, mask=rgba.getchannel("A"))
            tmp = tif + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(DPI, DPI), compression="tiff_lzw")
            os.replace(tmp, tif)
    png = os.path.join(PREV, os.path.basename(stem) + "_preview.png")
    fig.savefig(png, dpi=preview_dpi, bbox_inches=None, facecolor="white")
    print("saved", svg)
    print("saved", tif)
    print("preview", png)


def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    ct = A.obs["vendor_cluster"].astype(str).map(GROUPS).values

    counts = {g: int((ct == g).sum()) for g in GROUPS.values()}
    order = sorted(counts, key=lambda g: counts[g], reverse=True)
    print("true counts:", counts)
    print("display order:", order)

    P = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv"))
    P = P[P["rank"] <= 5].copy()
    P["grp_order"] = P["grp"].map({g: i for i, g in enumerate(order)})
    P = P.sort_values(["grp_order", "rank"]).drop(columns=["grp_order"]).reset_index(drop=True)

    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
    gi = [gpos[g] for g in P["gid"].astype(str)]
    E = A[:, gi].X.toarray().astype(np.float32)
    vmax = np.empty(E.shape[1], dtype=np.float32)
    for j in range(E.shape[1]):
        nz = E[:, j][E[:, j] > 0]
        vmax[j] = float(np.percentile(nz, 99)) if nz.size else 1.0
    Zall = np.clip(E / (vmax[None, :] + 1e-9), 0, 1)

    main_idx, marks = [], []
    for g in order:
        idx = np.where(ct == g)[0]
        main_idx.extend(idx.tolist())
        marks.append(len(main_idx))
    main_idx = np.asarray(main_idx, dtype=int)
    Z1 = Zall[main_idx].T

    small = ["MyoEpi", "M"]
    small_idx = np.concatenate([np.where(ct == g)[0] for g in small])
    Z2 = Zall[small_idx].T

    fig = plt.figure(figsize=(7.20, 5.00), facecolor="white")
    gs = GridSpec(3, 5, figure=fig,
                  width_ratios=[1.92, 3.88, 0.22, 1.05, 0.13],
                  height_ratios=[0.15, 0.13, 4.0],
                  left=0.012, right=0.985, top=0.955, bottom=0.135,
                  hspace=0.015, wspace=0.02)
    axl = fig.add_subplot(gs[2, 0])
    axb1, axn1, axh1 = [fig.add_subplot(gs[r, 1]) for r in range(3)]
    axb2, axn2, axh2 = [fig.add_subplot(gs[r, 3]) for r in range(3)]
    axc = fig.add_subplot(gs[2, 4])

    def draw_bar(axb, axn, nC, kind):
        bars = np.zeros((1, nC))
        if kind == "main":
            for i in range(len(order)):
                a = 0 if i == 0 else marks[i - 1]
                bars[0, a:marks[i]] = i
            cmap = ListedColormap([COL[g] for g in order])
            norm = BoundaryNorm(np.arange(-0.5, len(order) + 0.5), len(order))
        else:
            split = int((ct == "MyoEpi").sum())
            bars[0, :split] = 0
            bars[0, split:] = 1
            cmap = ListedColormap([COL["MyoEpi"], COL["M"]])
            norm = BoundaryNorm(np.arange(-0.5, 2.5), 2)
        axb.imshow(bars, aspect="auto", cmap=cmap, norm=norm,
                   extent=(0, nC, 0, 1), interpolation="nearest")
        axb.set_xticks([]); axb.set_yticks([])
        for sp in axb.spines.values():
            sp.set_visible(False)
        axn.set_xlim(0, nC); axn.set_ylim(0, 1); axn.axis("off")
        if kind == "small":
            split = int((ct == "MyoEpi").sum())
            axb.text(split / 2, 0.5, "MyoEpi", ha="center", va="center",
                     fontsize=5.0, weight="bold", color="white")
            axb.text(split + (nC - split) / 2, 0.5, "M", ha="center", va="center",
                     fontsize=5.0, weight="bold", color="white")
            axb.axvline(split, color="white", lw=0.7)
            axn.text(split / 2, 0.20, "n=65", ha="center", va="bottom", fontsize=5.0)
            axn.text(split + (nC - split) / 2, 0.20, "n=45", ha="center", va="bottom", fontsize=5.0)

    draw_bar(axb1, axn1, Z1.shape[1], "main")
    draw_bar(axb2, axn2, Z2.shape[1], "small")
    axb1.set_title("All groups  (n = 9,658 single cells)", fontsize=5.9, pad=3)
    axb2.set_title("MyoEpi + M (enlarged)", fontsize=5.9, pad=3)

    im = axh1.imshow(Z1, aspect="auto", cmap="viridis", vmin=0, vmax=1,
                     interpolation="antialiased", rasterized=True)
    axh2.imshow(Z2, aspect="auto", cmap="viridis", vmin=0, vmax=1,
                interpolation="nearest", rasterized=True)
    for axh in (axh1, axh2):
        axh.set_xticks([]); axh.set_yticks([])
        for sp in axh.spines.values():
            sp.set_visible(False)
    for i in range(1, len(order)):
        axh1.axvline(marks[i - 1], color="white", lw=0.55)
    split = int((ct == "MyoEpi").sum())
    axh2.axvline(split, color="white", lw=0.7)

    axl.set_xlim(0, 1)
    axl.set_ylim(Z1.shape[0] - 0.5, -0.5)
    axl.axis("off")
    for i, r in enumerate(P.itertuples()):
        gid = str(r.gid)
        short_id = gid[5:] if gid.startswith("ncbi_") else gid
        lab = str(r.sym) if str(r.sym) == short_id else "%s (%s)" % (r.sym, short_id)
        axl.text(1.0, i, lab, ha="right", va="center", fontsize=5.0)

    pos = axh2.get_position()
    cw, ch = 0.010, pos.height * 0.38
    axc.set_position([pos.x1 + 0.008, pos.y1 - ch - 0.010, cw, ch])
    cb = fig.colorbar(im, cax=axc)
    cb.set_label("Relative expression", fontsize=5.0, labelpad=2)
    cb.ax.tick_params(labelsize=5.0, length=1.0, pad=1.2)
    cb.outline.set_visible(False)

    # Matplotlib fills legend columns first; interleave to read row-wise.
    half = (len(order) + 1) // 2
    legend_order = []
    for i in range(half):
        legend_order.append(order[i])
        if i + half < len(order):
            legend_order.append(order[i + half])
    handles = [
        Patch(facecolor=COL[g], edgecolor="none",
              label="%s  n=%s" % (g, format(counts[g], ",")))
        for g in legend_order
    ]
    leg = fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.012),
        ncol=5, frameon=False, fontsize=5.0, handlelength=0.85,
        handletextpad=0.30, columnspacing=1.0, labelspacing=0.35,
        title=("All 9,658 cells; group width is proportional to the true cell number. "
               "Right: all MyoEpi and M cells."),
        title_fontsize=5.0,
    )
    leg.get_title().set_ha("center")

    stem = os.path.join(DELIV, "Fig2_singlecell_marker_heatmap_truecounts_highcontrast_inset")
    save_svg_tiff(fig, stem, preview_dpi=320)
    plt.close(fig)
    print("DONE", stem)


if __name__ == "__main__":
    log = open(os.path.join(ROOT, "logs", "fig2_truecounts_v3.log"), "w", encoding="utf8", buffering=1)
    sys.stdout = log
    sys.stderr = log
    main()
    log.close()

