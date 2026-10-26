# -*- coding: utf-8 -*-
"""213_fig5_marker_bubble_10groups.py
Marker-gene bubble plot for the ten current cell states.

Two interpretable, data-driven markers are selected per state
(20 genes total). Dot size = percent expressing; colour = scaled mean expression.
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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
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

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
ORDER = ["PH", "GH", "TS", "MH2", "F", "R", "B", "MyoEpi", "M", "PM-like"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "M": "#2166AC",
       "MyoEpi": "#92C5DE", "PM-like": "#7B3294"}
SELECTED = {
    "PH": ["PPAF2", "SSPO", "Grik2"],
    "GH": ["A2M", "CES2"],
    "TS": ["PEN-3", "APLP2"],
    "MH2": ["TIE1", "LOC113810024"],
    "F": ["Jon99Cii", "Chia"],
    "R": ["mesh", "Oacyl"],
    "B": ["MALRD1", "POLR2H"],
    "MyoEpi": ["lcc1", "KCNQ5"],
    "M": ["RASGRF1", "ADCY5"],
    "PM-like": ["Fas2", "LOC113809538"],
}


def flatten_tiff(path, dpi=1200):
    with Image.open(path) as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA")
            bg = Image.new("RGB", rgba.size, (255, 255, 255))
            bg.paste(rgba, mask=rgba.getchannel("A"))
            tmp = path + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(dpi, dpi), compression="tiff_lzw")
            os.replace(tmp, path)


def build_data():
    T = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv"))
    M = pd.read_csv(os.path.join(OUT, "markers_10groups.csv"))
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.obs["ct"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    var_names = set(A.var_names.astype(str))
    rows, selected_rows, gids = [], [], []
    for grp in ORDER:
        for gene in SELECTED[grp]:
            hit = T[(T.grp == grp) & (T.sym == gene)]
            if len(hit) != 1:
                raise ValueError(f"Expected one marker row for {grp} {gene}, got {len(hit)}")
            gid = str(hit.iloc[0].gid)
            if gid not in var_names:
                raise ValueError(f"Gene id not found in AnnData: {gid}")
            gids.append(gid)
            mh = M[(M.grp == grp) & (M.names.astype(str) == gid)]
            qval = float(mh["pvals_adj"].min()) if len(mh) else float("nan")
            pct_g = float(mh["pct_nz_group"].max()) if len(mh) else float("nan")
            pct_o = float(mh["pct_nz_reference"].max()) if len(mh) else float("nan")
            selected_rows.append({"group": grp, "gene": gene, "gid": gid,
                                  "marker_rank": int(hit.iloc[0]["rank"]),
                                  "log2FC": float(hit.iloc[0]["log2FC"]),
                                  "dPct": float(hit.iloc[0]["dPct"]),
                                  "pvals_adj": qval,
                                  "pct_group": pct_g,
                                  "pct_other": pct_o})
    E = A[:, gids].X.toarray().astype(np.float32)
    for j, (grp, gene) in enumerate([(g, x) for g in ORDER for x in SELECTED[g]]):
        x = E[:, j]
        for row_i, row_group in enumerate(ORDER):
            m = (A.obs["ct"].astype(str).values == row_group)
            rows.append({"group": row_group, "gene": gene, "gid": gids[j],
                         "pct_expressing": float(100 * (x[m] > 0).mean()),
                         "mean_logexpr": float(x[m].mean())})
    D = pd.DataFrame(rows)
    piv_mean = D.pivot(index="group", columns="gene", values="mean_logexpr").reindex(index=ORDER)
    piv_mean = piv_mean[[gene for grp in ORDER for gene in SELECTED[grp]]]
    z = (piv_mean - piv_mean.mean(axis=0)) / piv_mean.std(axis=0, ddof=0).replace(0, 1)
    D["zscore"] = [float(z.loc[r.group, r.gene]) for r in D.itertuples()]
    D.to_csv(os.path.join(OUT, "Fig5_marker_bubble_10groups_data.csv"),
             index=False, encoding="utf-8-sig")
    S = pd.DataFrame(selected_rows)
    S.to_csv(os.path.join(OUT, "Fig5_marker_bubble_10groups_selected.csv"),
             index=False, encoding="utf-8-sig")
    return D, S, z


def main():
    D, S, Z = build_data()
    genes = S["gene"].tolist()
    nG = len(genes)
    comp = pd.read_csv(os.path.join(OUT, "composition_new_names_10groups.csv"))
    counts = dict(zip(comp["group"].astype(str), comp["n_total"].astype(int)))
    group_pos = {}
    start = 0
    for grp in ORDER:
        stop = start + len(SELECTED[grp]) - 1
        group_pos[grp] = (start, stop, (start + stop) / 2.0)
        start = stop + 1
    separators = [group_pos[grp][1] + 0.5 for grp in ORDER[:-1]]
    fig = plt.figure(figsize=(6.40, 5.5 / 2.54), facecolor="white")
    ax = fig.add_axes([0.105, 0.290, 0.700, 0.570])
    sx = np.arange(nG)
    sy = np.arange(len(ORDER))
    pct = D.pivot(index="group", columns="gene", values="pct_expressing").reindex(index=ORDER)[genes]
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=-2, vmax=2)
    for i, grp in enumerate(ORDER):
        for j, gene in enumerate(genes):
            v = float(Z.loc[grp, gene])
            p = float(pct.loc[grp, gene])
            ax.scatter(j, i, s=5 + 62 * p / 100.0, c=[cmap(norm(v))],
                       edgecolors="white", linewidths=0.35, zorder=3)
    # subtle grid and block separators
    for i in range(len(ORDER)):
        ax.axhline(i, color="#E9E9E9", lw=0.35, zorder=0)
    for j in range(nG):
        ax.axvline(j, color="#F0F0F0", lw=0.30, zorder=0)
    for sep in separators:
        ax.axvline(sep, color="#BDBDBD", lw=0.55, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(-0.6, nG - 0.4)
    ax.set_ylim(len(ORDER) - 0.55, -0.55)
    ax.set_xticks(sx)
    ax.set_xticklabels(genes, rotation=48, ha="right", rotation_mode="anchor", fontsize=5.0)
    ax.set_yticks(sy)
    ax.set_yticklabels([f"{g}  (n={counts.get(g, 0):,})" for g in ORDER], fontsize=5.0)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", length=0, pad=3)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # Coloured strips mark gene blocks; row labels already name the cell states.
    for grp in ORDER:
        a, b, _ = group_pos[grp]
        ax.plot([a - 0.35, b + 0.35], [-0.42, -0.42], color=COL[grp],
                lw=2.2, solid_capstyle="butt", zorder=4)
    ax.set_title("Marker genes across the ten cell states", fontsize=6.5, pad=8)
    # size legend in a separate inset, below the colour bar
    ax_leg = fig.add_axes([0.835, 0.18, 0.135, 0.22])
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1); ax_leg.axis("off")
    ax_leg.text(0.5, 0.95, "% expressing", ha="center", va="top", fontsize=5.0,
                color="#333333")
    for p, y in zip([25, 50, 75], [0.68, 0.43, 0.18]):
        ax_leg.scatter([0.26], [y], s=5 + 62 * p / 100.0, c="white",
                       edgecolors="#555555", linewidths=0.5)
        ax_leg.text(0.52, y, f"{p}%", ha="left", va="center", fontsize=5.0,
                    color="#333333")
    # colour bar
    cax = fig.add_axes([0.875, 0.50, 0.018, 0.34])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("Scaled mean expression", fontsize=5.2, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.5, width=0.4)
    cb.outline.set_linewidth(0.4)
    stem = os.path.join(PARTS, "Fig5_marker_bubble_10groups")
    bbox = fig.bbox_inches
    fig.savefig(stem + ".svg", bbox_inches=bbox, pad_inches=0, facecolor="white")
    fig.savefig(stem + "_qa.pdf", bbox_inches=bbox, pad_inches=0, facecolor="white")
    fig.savefig(stem + ".tiff", dpi=1200, bbox_inches=bbox, pad_inches=0,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    flatten_tiff(stem + ".tiff", 1200)
    fig.savefig(stem + "_preview.png", dpi=450, bbox_inches=bbox, pad_inches=0,
                facecolor="white")
    print("saved", stem + ".svg")
    print("saved", stem + ".tiff")
    print("saved", stem + "_preview.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
