# -*- coding: utf-8 -*-
"""222_ph_dpt_heatmap_preview.py
Pseudotime-ordered marker heatmap for PH substates using DPT.
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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
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
COL = {"PH_s0": "#59A14F", "PH_s1": "#E15759", "PH_s2": "#F28E2B", "PH_s3": "#4E79A7"}
GENES = ["PPAF2", "Igfbp7", "SSPO", "Grik2", "ccnd2", "TGM1", "SERPINB1", "slo",
         "Syne1", "LOC113808957", "hdc", "LOC113810024", "PEN-3", "THBS4"]


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
    ph = A[A.obs["ct"].astype(str) == "PH"].copy()
    ph.obs["sub"] = sub.reindex(ph.obs_names)["sub"].values
    ph.obs["group2"] = ph.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    ph.X = ph.layers["counts"].copy()
    sc.pp.normalize_total(ph, target_sum=1e4)
    sc.pp.log1p(ph)

    M = pd.read_csv(os.path.join(OUT, "markers_10groups.csv"))
    M = M[M.grp == "PH"].sort_values("scores", ascending=False)
    sym2gid = {}
    for r in M.itertuples():
        sym2gid.setdefault(r.sym, r.names)
    gids = []
    for g in GENES:
        if g not in sym2gid:
            raise ValueError(f"Missing PH gene {g}")
        gids.append(sym2gid[g])
    E = ph[:, gids].X.toarray().astype(float)
    pt = ph.obs["dpt_pseudotime"].astype(float).values
    bins = pd.qcut(pt, 20, labels=False, duplicates="drop")
    nb = int(np.nanmax(bins)) + 1
    Z = np.zeros((len(GENES), nb))
    raw = np.zeros_like(Z)
    dpt_med = np.zeros(nb)
    dominant = []
    for b in range(nb):
        m = bins == b
        raw[:, b] = E[m].mean(0)
        dpt_med[b] = np.median(pt[m])
        vc = ph.obs.loc[m, "sub"].astype(str).value_counts()
        dominant.append(vc.index[0] if len(vc) else "PH_s0")
    for j in range(len(GENES)):
        mu, sd = raw[j].mean(), raw[j].std(ddof=0)
        Z[j] = (raw[j] - mu) / sd if sd > 1e-9 else 0
    # save data
    rows = []
    for j, g in enumerate(GENES):
        for b in range(nb):
            rows.append({"gene": g, "bin": b + 1, "median_dpt": dpt_med[b],
                         "mean_logexpr": raw[j, b], "zscore": Z[j, b],
                         "dominant_substate": dominant[b]})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "Fig5_PH_dpt_heatmap_data.csv"),
                              index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.0, 4.2), facecolor="white")
    ax_top = fig.add_axes([0.22, 0.80, 0.64, 0.045])
    ax = fig.add_axes([0.22, 0.17, 0.64, 0.61])
    # top dominant-substate strip
    idx = [PH_KEYS.index("PH_s" + d.split("_s")[1]) if d.startswith("PH_s") else 0 for d in dominant]
    ax_top.imshow(np.array([idx]), aspect="auto", cmap=ListedColormap([COL[k] for k in PH_KEYS]),
                  vmin=-0.5, vmax=len(PH_KEYS) - 0.5)
    ax_top.set_xticks([]); ax_top.set_yticks([])
    for sp in ax_top.spines.values():
        sp.set_visible(False)
    im = ax.imshow(np.clip(Z, -2, 2), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(np.arange(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=5.2)
    ax.tick_params(axis="y", length=0, pad=2)
    xt = [0, 5, 10, 15, nb - 1]
    ax.set_xticks(xt)
    ax.set_xticklabels([str(k + 1) for k in xt], fontsize=5.0)
    ax.set_xlabel("DPT pseudotime bins (early → late)", fontsize=5.5)
    ax.tick_params(axis="x", length=1.3, width=0.4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = fig.add_axes([0.875, 0.26, 0.018, 0.42])
    cb = fig.colorbar(im, cax=cax); cb.set_label("Scaled mean expression", fontsize=5.2, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.3, width=0.4); cb.outline.set_linewidth(0.4)
    handles = [Patch(color=COL[k], label=k.split("_")[1]) for k in PH_KEYS]
    ax_top.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.05), ncol=4,
                  frameon=False, fontsize=5.0, handlelength=0.7, columnspacing=0.8,
                  handletextpad=0.25)
    ax_top.set_title("Dominant PH substate along pseudotime", fontsize=5.6, pad=14)
    fig.text(0.50, 0.955, "PH pseudotime-ordered marker dynamics", ha="center", va="top", fontsize=7.2)
    fig.text(0.50, 0.025, "20 equal-count DPT bins; colour = row z-scored mean log expression. "
             "Bins are ordered by DPT, not equal time intervals.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    save_export(fig, "Fig5_PH_DPT_pseudotime_heatmap")
    plt.close(fig)


if __name__ == "__main__":
    main()
