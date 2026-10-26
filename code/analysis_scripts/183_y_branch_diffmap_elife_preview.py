# -*- coding: utf-8 -*-
"""183_y_branch_diffmap_elife_preview.py
Reference-style branched hemocyte map: PH -> TS -> GH / MH2.
Balanced diffusion components; TS downsampled only for visualization.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse, FancyArrowPatch
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from PIL import Image

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
MAP = {"MH1": "TS", "PH": "PH", "GH": "GH", "MH2": "MH2"}
ORDER = ["PH", "GH", "TS", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A"}
LIN_LAB = {"PH": "Prohemocyte", "GH": "Granulocyte lineage", "MH2": "Monocytic hemocyte lineage"}

def ell(pts, nsd=1.65):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]
    ev, evec = ev[o], evec[:, o]
    ang = np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), ang

def bezier(p0, p2, bend=0.15, n=160):
    p0 = np.asarray(p0, float); p2 = np.asarray(p2, float)
    mid = (p0 + p2) / 2.0
    v = p2 - p0
    perp = np.array([-v[1], v[0]])
    p1 = mid + perp * bend
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t ** 2 * p2


def draw_arrowhead(ax, curve, color, lw=1.8):
    p_end = curve[-1]; p_prev = curve[-8]
    ax.add_patch(FancyArrowPatch(p_prev, p_end, arrowstyle="-|>,head_length=0.55,head_width=0.26",
                                 mutation_scale=11, lw=lw, color=color, alpha=0.98, zorder=9,
                                 shrinkA=0, shrinkB=0,
                                 path_effects=[pe.withStroke(linewidth=2.8, foreground="white")]))


def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).map(MAP).values
    rng = np.random.default_rng(7)
    ts = np.where(ct == "TS")[0]
    keep = np.where(ct != "TS")[0].tolist() + rng.choice(ts, size=min(len(ts), 900), replace=False).tolist()
    B = A[np.array(sorted(keep))].copy()
    B.obs["ct_new"] = pd.Categorical(B.obs["ct"].astype(str).map(MAP), categories=ORDER)
    sc.pp.neighbors(B, use_rep="X_pca_harmony", n_neighbors=20, n_pcs=30)
    sc.tl.diffmap(B, n_comps=15)
    D = np.asarray(B.obsm["X_diffmap"])[:, 1:3].copy()
    D[:, 0] *= -1  # orient PH lower-left, GH right
    bct = B.obs["ct_new"].astype(str).values

    fig = plt.figure(figsize=(4.0, 3.7), facecolor="white")
    ax = fig.add_axes([0.12, 0.13, 0.76, 0.74])
    for g in ORDER:
        m = bct == g
        ax.scatter(D[m, 0], D[m, 1], s=1.3, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%s)" % (g, format(int(m.sum()), ",")))
    cent = {g: D[bct == g].mean(0) for g in ORDER}
    # smooth shared trunk and branches
    trunk = bezier(cent["PH"], cent["TS"], bend=0.06)
    br_gh = bezier(cent["TS"], cent["GH"], bend=0.20)
    br_mh2 = bezier(cent["TS"], cent["MH2"], bend=-0.18)
    for curve, color, lw in [(trunk, "#333333", 2.0),
                             (br_gh, COL["GH"], 2.0),
                             (br_mh2, COL["MH2"], 2.0)]:
        ax.plot(curve[:, 0], curve[:, 1], lw=lw, color=color, zorder=8,
                solid_capstyle="round", path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])
        draw_arrowhead(ax, curve, color, lw=lw)
    # labels
    ax.annotate(LIN_LAB["PH"], xy=cent["PH"], xytext=(-16, -16), textcoords="offset points",
                fontsize=5.5, color=COL["PH"], weight="bold", ha="center", va="top",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=9)
    ax.annotate(LIN_LAB["GH"], xy=cent["GH"], xytext=(20, -12), textcoords="offset points",
                fontsize=5.5, color=COL["GH"], weight="bold", ha="center", va="top",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=9)
    ax.annotate(LIN_LAB["MH2"], xy=cent["MH2"], xytext=(6, 16), textcoords="offset points",
                fontsize=5.5, color=COL["MH2"], weight="bold", ha="left", va="bottom",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=9)
    ax.annotate("TS", xy=cent["TS"], xytext=(0, -13), textcoords="offset points",
                fontsize=5.4, color="#333333", weight="bold", ha="center", va="top",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("Component 1", fontsize=6.4); ax.set_ylabel("Component 2", fontsize=6.4)
    for sp in ax.spines.values():
        sp.set_color("#777777"); sp.set_linewidth(0.45)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=4, frameon=False,
              markerscale=4.2, fontsize=5.0, handletextpad=0.25, columnspacing=0.8)
    fig.text(0.50, 0.018,
             "Balanced diffusion map; TS shown as 900 cells for visualization.",
             ha="center", va="bottom", fontsize=5.0, color="#333333")
    p = os.path.join(PREV, "02d_branch_diffmap_elife_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    # final exports
    deliv = os.path.join(r"D:\桌面\single cell\文章配图\最终_2026-10-17", "拟时分叉图")
    os.makedirs(deliv, exist_ok=True)
    stem = os.path.join(deliv, "Fig_pseudotime_branch_diffmap")
    fig.savefig(stem + ".svg", bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(stem + ".tiff", dpi=1200, bbox_inches="tight", pad_inches=0.04,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(stem + ".tiff") as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA"); bg = Image.new("RGB", rgba.size, (255,255,255))
            bg.paste(rgba, mask=rgba.getchannel("A")); tmp = stem + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(1200,1200), compression="tiff_lzw")
            os.replace(tmp, stem + ".tiff")
    print("saved", stem + ".svg")
    print("saved", stem + ".tiff")
    plt.close(fig)

if __name__ == "__main__":
    main()
