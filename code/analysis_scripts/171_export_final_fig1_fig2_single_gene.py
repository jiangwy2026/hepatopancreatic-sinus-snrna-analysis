# -*- coding: utf-8 -*-
"""171_export_final_fig1_fig2_single_gene.py
Final exports requested on 2026-10-17:
- keep the original combined t-SNE;
- add separately rendered control and hemospasia t-SNE panels;
- Figure 2 as a single-cell marker heatmap with an enlarged M + MyoEpi block;
- save every single-gene t-SNE as its own editable SVG and 1200-dpi TIFF.
Only SVG and TIFF are placed in the delivery folder.
"""
import os, re, math, json
import matplotlib
import matplotlib as mpl
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image
import sys

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
ORDER = ["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]
FULL = {"PH": "Prohemocyte", "GH": "Granular hemocyte", "TS": "Biosynthetic transition",
        "MH2": "Macrophage-like", "F": "F cell", "R": "R cell", "B": "B cell",
        "M": "M cell", "MyoEpi": "Myoepithelial-like", "PM-like": "Peritrophic matrix-like"}
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "M": "#2166AC",
       "MyoEpi": "#92C5DE", "PM-like": "#7B3294"}
SAMPLE_TITLE = {"blank": "Control (0 h)", "treated": "Hemospasia (6 h)"}
GREY = "#D0D0D0"
GCMAP = ListedColormap(plt.cm.Greens(np.linspace(0.25, 0.95, 256)))
DPI = 1200


def save_svg_tiff(fig, stem, preview=True, preview_dpi=300):
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
    if preview:
        os.makedirs(PREV, exist_ok=True)
        png = os.path.join(PREV, os.path.basename(stem) + "_preview.png")
        fig.savefig(png, dpi=preview_dpi, bbox_inches=None, facecolor="white")
    print("saved", svg)
    print("saved", tif)


def load_data():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))[T.obs_names].copy()
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    ct = A.obs["vendor_cluster"].astype(str).map(GROUPS).values
    sample = A.obs["sample"].astype(str).values
    return T, A, tsne, ct, sample


def common_limits(tsne, pad_frac=0.035):
    xmin, xmax = float(np.min(tsne[:, 0])), float(np.max(tsne[:, 0]))
    ymin, ymax = float(np.min(tsne[:, 1])), float(np.max(tsne[:, 1]))
    px = (xmax - xmin) * pad_frac
    py = (ymax - ymin) * pad_frac
    return (xmin - px, xmax + px), (ymin - py, ymax + py)


def set_tsne_axes(ax, xlim, ylim, axis_labels=True):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_box_aspect(1)
    if axis_labels:
        ax.set_xlabel("t-SNE 1", fontsize=6.4, labelpad=2)
        ax.set_ylabel("t-SNE 2", fontsize=6.4, labelpad=2)
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
    for sp in ax.spines.values():
        sp.set_color("#8F8F8F")
        sp.set_linewidth(0.45)


def draw_tsne(ax, tsne, ct, xlim, ylim, mask=None, label_min=30):
    if mask is None:
        mask = np.ones(tsne.shape[0], dtype=bool)
    for g in ORDER:
        m = mask & (ct == g)
        n = int(m.sum())
        if n:
            ax.scatter(tsne[m, 0], tsne[m, 1], s=1.45, c=COL[g], linewidths=0,
                       rasterized=True, label="%s (%s)" % (g, format(n, ",")))
    for g in ORDER:
        m = mask & (ct == g)
        if int(m.sum()) < label_min:
            continue
        ax.text(float(np.median(tsne[m, 0])), float(np.median(tsne[m, 1])), g,
                fontsize=5.7, weight="bold", ha="center", va="center", color="#111111",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=8)
    set_tsne_axes(ax, xlim, ylim)


def fig1_original(T, tsne, ct):
    xlim, ylim = common_limits(tsne)
    fig = plt.figure(figsize=(5.65, 4.25), facecolor="white")
    ax = fig.add_axes([0.105, 0.13, 0.565, 0.79])
    for g in ORDER:
        m = ct == g
        ax.scatter(tsne[m, 0], tsne[m, 1], s=1.55, c=COL[g], linewidths=0,
                   rasterized=True, label="%s  %s (%s)" % (g, FULL[g], format(int(m.sum()), ",")))
    for g in ORDER:
        m = ct == g
        ax.text(float(np.median(tsne[m, 0])), float(np.median(tsne[m, 1])), g,
                fontsize=5.7, weight="bold", ha="center", va="center", color="#111111",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=8)
    set_tsne_axes(ax, xlim, ylim)
    ax.legend(loc="center left", bbox_to_anchor=(1.015, 0.5), frameon=False,
              markerscale=4.4, fontsize=5.45, handletextpad=0.30, labelspacing=0.32)
    ax.set_title("All %s cells" % format(T.n_obs, ","), fontsize=6.7, pad=5)
    save_svg_tiff(fig, os.path.join(DELIV, "Fig1b_TSNE_all_groups"), preview=True)
    plt.close(fig)


def fig1_by_group(tsne, ct, sample):
    xlim, ylim = common_limits(tsne)
    for smp, short in (("blank", "control"), ("treated", "hemospasia")):
        m0 = sample == smp
        fig = plt.figure(figsize=(3.50, 4.05), facecolor="white")
        ax = fig.add_axes([0.13, 0.285, 0.74, 0.64])
        draw_tsne(ax, tsne, ct, xlim, ylim, mask=m0, label_min=30)
        ax.set_title("%s   n = %s" % (SAMPLE_TITLE[smp], format(int(m0.sum()), ",")),
                     fontsize=6.8, pad=6)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.155), ncol=2,
                  frameon=False, markerscale=4.2, fontsize=5.15,
                  handletextpad=0.25, columnspacing=0.75, labelspacing=0.28)
        stem = os.path.join(DELIV, "Fig1c_TSNE_%s" % short)
        save_svg_tiff(fig, stem, preview=True)
        plt.close(fig)


def fig1_comparison(tsne, ct, sample):
    xlim, ylim = common_limits(tsne)
    fig = plt.figure(figsize=(7.20, 3.90), facecolor="white")
    axes = [fig.add_axes([0.055, 0.14, 0.37, 0.74]),
            fig.add_axes([0.455, 0.14, 0.37, 0.74])]
    for ax, smp in zip(axes, ("blank", "treated")):
        m0 = sample == smp
        draw_tsne(ax, tsne, ct, xlim, ylim, mask=m0, label_min=30)
        ax.set_title("%s   n = %s" % (SAMPLE_TITLE[smp], format(int(m0.sum()), ",")),
                     fontsize=6.7, pad=6)
    handles, labels = axes[0].get_legend_handles_labels()
    n_blank = {g: int(((sample == "blank") & (ct == g)).sum()) for g in ORDER}
    n_treated = {g: int(((sample == "treated") & (ct == g)).sum()) for g in ORDER}
    labels = ["%s (%s / %s)" % (g, format(n_blank[g], ","), format(n_treated[g], ",")) for g in ORDER]
    axes[1].legend(handles, labels, loc="center left", bbox_to_anchor=(1.045, 0.5),
                   frameon=False, markerscale=4.3, fontsize=5.35,
                   handletextpad=0.30, labelspacing=0.32,
                   title="Control / Hemospasia", title_fontsize=5.35)
    save_svg_tiff(fig, os.path.join(DELIV, "Fig1c_TSNE_control_and_hemospasia"), preview=True)
    plt.close(fig)


def fig2_heatmap(A, ct):
    P = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv"))
    P = P[P["rank"] <= 5].copy().reset_index(drop=True)
    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
    gi = [gpos[g] for g in P["gid"].astype(str)]
    E = A[:, gi].X.toarray().astype(np.float32)
    vmax = np.empty(E.shape[1], dtype=np.float32)
    for j in range(E.shape[1]):
        nz = E[:, j][E[:, j] > 0]
        vmax[j] = float(np.percentile(nz, 99)) if nz.size else 1.0
    Zall = np.clip(E / (vmax[None, :] + 1e-9), 0, 1)

    rng = np.random.default_rng(0)
    main_idx, marks, shown, total = [], [], [], []
    for g in ORDER:
        idx = np.where(ct == g)[0]
        n = len(idx)
        k = min(n, max(15, int(round(25 * math.log10(max(1, n))))))
        pick = np.sort(rng.choice(idx, k, replace=False))
        main_idx.extend(pick.tolist())
        marks.append(len(main_idx))
        shown.append(k)
        total.append(n)
    main_idx = np.asarray(main_idx, dtype=int)
    small = ["M", "MyoEpi"]
    small_idx = np.concatenate([np.where(ct == g)[0] for g in small])
    Z1 = Zall[main_idx].T
    Z2 = Zall[small_idx].T

    fig = plt.figure(figsize=(7.20, 5.15), facecolor="white")
    gs = GridSpec(3, 5, figure=fig,
                  width_ratios=[2.00, 3.45, 0.28, 1.52, 0.17],
                  height_ratios=[0.15, 0.13, 4.0],
                  left=0.012, right=0.985, top=0.955, bottom=0.065,
                  hspace=0.015, wspace=0.02)
    axl = fig.add_subplot(gs[2, 0])
    axb1, axn1, axh1 = [fig.add_subplot(gs[r, 1]) for r in range(3)]
    axb2, axn2, axh2 = [fig.add_subplot(gs[r, 3]) for r in range(3)]
    axc = fig.add_subplot(gs[2, 4])

    def draw_panel(axb, axn, axh, Z, marks_, total_, kind):
        nC = Z.shape[1]
        if kind == "main":
            bars = np.zeros((1, nC))
            for i in range(len(ORDER)):
                a = 0 if i == 0 else marks_[i - 1]
                bars[0, a:marks_[i]] = i
            cmap = ListedColormap([COL[g] for g in ORDER])
            norm = BoundaryNorm(np.arange(-0.5, len(ORDER) + 0.5), len(ORDER))
        else:
            a = int((ct == "M").sum())
            bars = np.zeros((1, nC))
            bars[0, :a] = 0
            bars[0, a:] = 1
            cmap = ListedColormap([COL["M"], COL["MyoEpi"]])
            norm = BoundaryNorm(np.arange(-0.5, 2.5), 2)
        axb.imshow(bars, aspect="auto", cmap=cmap, norm=norm,
                   extent=(0, nC, 0, 1), interpolation="nearest")
        axb.set_xticks([]); axb.set_yticks([])
        for sp in axb.spines.values():
            sp.set_visible(False)
        im = axh.imshow(Z, aspect="auto", cmap="magma", vmin=0, vmax=1,
                        interpolation="nearest", rasterized=True)
        axh.set_xticks([]); axh.set_yticks([])
        for sp in axh.spines.values():
            sp.set_visible(False)
        axn.set_xlim(0, nC); axn.set_ylim(0, 1); axn.axis("off")
        if kind == "main":
            for i, g in enumerate(ORDER):
                a = 0 if i == 0 else marks_[i - 1]
                b = marks_[i]
                axb.text((a + b) / 2, 0.5, g, ha="center", va="center",
                         fontsize=5.0, weight="bold", color="white")
                if i > 0:
                    axb.axvline(a, color="white", lw=0.65)
                    axh.axvline(a, color="white", lw=0.65)
                axn.text((a + b) / 2, 0.20, "n=%s" % format(total_[i], ","),
                         ha="center", va="bottom", fontsize=5.0)
        else:
            a = int((ct == "M").sum())
            axb.text(a / 2, 0.5, "M", ha="center", va="center",
                     fontsize=5.0, weight="bold", color="white")
            axb.text(a + (nC - a) / 2, 0.5, "MyoEpi", ha="center", va="center",
                     fontsize=5.0, weight="bold", color="white")
            axb.axvline(a, color="white", lw=0.65)
            axh.axvline(a, color="white", lw=0.65)
            axn.text(a / 2, 0.20, "n=%s" % format(int((ct == "M").sum()), ","),
                     ha="center", va="bottom", fontsize=5.0)
            axn.text(a + (nC - a) / 2, 0.20,
                     "n=%s" % format(int((ct == "MyoEpi").sum()), ","),
                     ha="center", va="bottom", fontsize=5.0)
        return im

    im = draw_panel(axb1, axn1, axh1, Z1, marks, total, "main")
    draw_panel(axb2, axn2, axh2, Z2, None, None, "small")
    axb1.set_title("All groups", fontsize=5.9, pad=3)
    axb2.set_title("M + MyoEpi (enlarged)", fontsize=5.9, pad=3)

    axl.set_xlim(0, 1)
    axl.set_ylim(Z1.shape[0] - 0.5, -0.5)
    axl.axis("off")
    for i, r in enumerate(P.itertuples()):
        gid = str(r.gid)
        short_id = gid[5:] if gid.startswith("ncbi_") else gid
        if str(r.sym) == short_id:
            lab = str(r.sym)
        else:
            lab = "%s (%s)" % (r.sym, short_id)
        axl.text(1.0, i, lab, ha="right", va="center", fontsize=5.0)

    pos = axh2.get_position()
    cw, ch = 0.012, pos.height * 0.38
    axc.set_position([pos.x1 + 0.010, pos.y1 - ch - 0.010, cw, ch])
    cb = fig.colorbar(im, cax=axc)
    cb.set_label("Relative expression", fontsize=5.0, labelpad=2)
    cb.ax.tick_params(labelsize=5.0, length=1.0, pad=1.2)
    cb.outline.set_visible(False)
    fig.text(0.50, 0.018,
             "Main: individual cells; group width scaled by log10(n). Inset: all M and MyoEpi cells.",
             ha="center", va="bottom", fontsize=5.0, color="#3A3A3A")
    save_svg_tiff(fig, os.path.join(DELIV, "Fig2_singlecell_marker_heatmap_small_cluster_inset"),
                  preview=True, preview_dpi=320)
    plt.close(fig)


def single_gene_tiles(A, tsne, ct):
    picks = pd.read_csv(os.path.join(OUT, "single_gene_panel_10groups.csv"))
    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
    xlim, ylim = common_limits(tsne)
    gene_dir = os.path.join(DELIV, "FigS_single_gene_TSNE")
    os.makedirs(gene_dir, exist_ok=True)
    rows = []
    for k, r in enumerate(picks.itertuples(), start=1):
        gid = str(r.gid)
        sym = str(r.gene)
        grp = str(r.group)
        v = A[:, gpos[gid]].X.toarray().ravel().astype(np.float32)
        pos = np.where(v > 0)[0]
        pct_group = 100.0 * float((v[ct == grp] > 0).mean())
        pct_all = 100.0 * float((v > 0).mean())
        fig = plt.figure(figsize=(2.70, 2.95), facecolor="white")
        ax = fig.add_axes([0.145, 0.175, 0.665, 0.61])
        ax.scatter(tsne[:, 0], tsne[:, 1], c=GREY, s=0.72, linewidths=0,
                   rasterized=True, zorder=1)
        if pos.size:
            vmax = float(np.percentile(v[pos], 95))
            rel = np.clip(v[pos] / (vmax + 1e-9), 0, 1)
            o = np.argsort(rel)
            sm = ax.scatter(tsne[pos[o], 0], tsne[pos[o], 1], c=rel[o],
                            cmap=GCMAP, vmin=0, vmax=1, s=1.55,
                            linewidths=0, rasterized=True, zorder=2)
            cax = fig.add_axes([0.845, 0.315, 0.030, 0.33])
            cb = fig.colorbar(sm, cax=cax)
            cb.set_ticks([0, 0.5, 1.0])
            cb.ax.tick_params(labelsize=5.0, length=0.9, pad=1.0)
            cb.outline.set_visible(False)
            cb.set_label("Relative expression", fontsize=5.0, labelpad=2)
        set_tsne_axes(ax, xlim, ylim, axis_labels=True)
        ax.set_title("%s  |  %s\n%.1f%% of group; %.1f%% of all cells" %
                     (sym, grp, pct_group, pct_all),
                     fontsize=5.55, pad=4, linespacing=1.12)
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", sym)
        stem = os.path.join(gene_dir, "%02d_%s_%s" % (k, grp, safe))
        save_svg_tiff(fig, stem, preview=True, preview_dpi=300)
        plt.close(fig)
        rows.append({"index": k, "group": grp, "gene": sym, "gid": gid,
                     "pct_group": round(pct_group, 3), "pct_all": round(pct_all, 3),
                     "svg": os.path.basename(stem) + ".svg",
                     "tiff": os.path.basename(stem) + ".tiff"})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "final_single_gene_export_manifest.csv"),
                              index=False, encoding="utf-8-sig")
    print("individual gene files:", len(rows))


def main():
    T, A, tsne, ct, sample = load_data()
    print("cells:", T.n_obs, "blank:", int((sample == "blank").sum()),
          "treated:", int((sample == "treated").sum()))
    fig1_original(T, tsne, ct)
    fig1_by_group(tsne, ct, sample)
    fig1_comparison(tsne, ct, sample)
    fig2_heatmap(A, ct)
    single_gene_tiles(A, tsne, ct)
    print("DONE", DELIV)


if __name__ == "__main__":
    log = open(os.path.join(ROOT, "logs", "final_export_171.log"), "w", encoding="utf8", buffering=1)
    sys.stdout = log
    sys.stderr = log
    main()
    log.close()


