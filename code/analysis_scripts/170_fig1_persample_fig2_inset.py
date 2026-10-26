# -*- coding: utf-8 -*-
"""170_fig1_persample_fig2_inset.py
(1) Fig1 增加 control / hemospasia 两张独立 t-SNE（原 UMAP、合并 t-SNE 保留）
(2) Fig2 单细胞热图 + 右侧小群（M, MyoEpi）放大块
"""
import io, sys, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
GROUPS = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER = ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
FULL = {"PH":"Prohemocyte","GH":"Granular hemocyte","TS":"Biosynthetic transition","MH2":"Macrophage-like",
        "F":"F cell","R":"R cell","B":"B cell","M":"M cell","MyoEpi":"Myoepithelial-like","PM-like":"Peritrophic matrix-like"}
COL = {"PH":"#B2182B","GH":"#E08214","TS":"#D6604D","MH2":"#8C510A","F":"#006D2C","R":"#41AB5D",
       "B":"#78C679","M":"#2166AC","MyoEpi":"#92C5DE","PM-like":"#7B3294"}

def fig1_persample():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    ct = T.obs["vendor_cluster"].astype(str).map(GROUPS).values
    s = T.obs["sample"].astype(str).values
    U = np.asarray(T.obsm["X_tsne"])
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.9))
    for ax, smp, title in zip(axes, ("blank", "treated"), ("Control (0 h)", "Hemospasia (6 h)")):
        m0 = s == smp
        for g in ORDER:
            m = m0 & (ct == g)
            if m.sum(): ax.scatter(U[m, 0], U[m, 1], s=1.4, c=COL[g], linewidths=0, rasterized=True,
                                   label="%s (%d)" % (g, int(m.sum())))
        for g in ORDER:
            m = m0 & (ct == g)
            if m.sum() < 30: continue
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), g, fontsize=5.8, weight="bold", ha="center",
                    va="center", color="black", path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=6)
        ax.set_title("%s   n = %s" % (title, format(int(m0.sum()), ",")), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
                   fontsize=5.8, handletextpad=0.3, labelspacing=0.28)
    p = os.path.join(FIG, "Fig1c_TSNE_control_vs_hemospasia_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)
    # 两张单独文件
    for smp, name in (("blank", "control"), ("treated", "hemospasia")):
        fig, ax = plt.subplots(figsize=(3.9, 3.7))
        m0 = s == smp
        for g in ORDER:
            m = m0 & (ct == g)
            if m.sum(): ax.scatter(U[m, 0], U[m, 1], s=1.4, c=COL[g], linewidths=0, rasterized=True,
                                   label="%s (%d)" % (g, int(m.sum())))
        for g in ORDER:
            m = m0 & (ct == g)
            if m.sum() < 30: continue
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), g, fontsize=5.8, weight="bold", ha="center",
                    va="center", color="black", path_effects=[pe.withStroke(linewidth=1.7, foreground="white")], zorder=6)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False, markerscale=4,
                  fontsize=5.6, handletextpad=0.25, columnspacing=0.8)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
        p = os.path.join(FIG, "Fig1c_TSNE_%s_preview.png" % name)
        fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

def fig2_inset():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.obs["ct"] = pd.Categorical(A.obs["vendor_cluster"].astype(str).map(GROUPS), categories=ORDER)
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    P = pd.read_csv(os.path.join(OUT, "top_markers_10groups.csv")); P = P[P["rank"] <= 5].copy()
    gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
    rng = np.random.default_rng(0)
    small = ["M", "MyoEpi"]
    cell_idx, marks, shown, total = [], [], [], []
    for g in ORDER:
        idx = np.where((A.obs["ct"] == g).values)[0]
        n = len(idx); k = min(n, max(15, int(round(25 * math.log10(n)))))
        pick = np.sort(rng.choice(idx, k, replace=False))
        cell_idx += list(pick); marks.append(len(cell_idx)); shown.append(k); total.append(n)
    main_idx = np.array(cell_idx)
    small_idx = np.concatenate([np.where((A.obs["ct"] == g).values)[0] for g in small])
    gi = [gpos[g] for g in P["gid"].astype(str)]
    def zmat(idx):
        E = A[idx][:, gi].X.toarray().astype(np.float32)
        return np.clip(E / (np.percentile(E, 99, axis=0) + 1e-9), 0, 1).T
    Z1, Z2 = zmat(main_idx), zmat(small_idx)

    fig = plt.figure(figsize=(9.8, 5.3))
    gs = GridSpec(3, 5, height_ratios=[0.030, 0.020, 1.0], width_ratios=[0.24, 1.0, 0.10, 0.34, 0.024],
                  hspace=0.02, wspace=0.04, left=0.02, right=0.96, top=0.965, bottom=0.03)
    axb1 = fig.add_subplot(gs[0, 1]); axn1 = fig.add_subplot(gs[1, 1]); axh1 = fig.add_subplot(gs[2, 1])
    axb2 = fig.add_subplot(gs[0, 3]); axn2 = fig.add_subplot(gs[1, 3]); axh2 = fig.add_subplot(gs[2, 3])
    axl = fig.add_subplot(gs[2, 0]); axc = fig.add_subplot(gs[2, 4])
    for axb, axn, axh, Z, marks_, total_, label in (
            (axb1, axn1, axh1, Z1, marks, total, "all groups"),
            (axb2, axn2, axh2, Z2, [len(Z2[0])], [len(small_idx)], "M + MyoEpi (enlarged)")):
        nC = Z.shape[1]
        bars = np.zeros((1, nC))
        if axb is axb1:
            for i in range(len(ORDER)):
                a = 0 if i == 0 else marks_[i - 1]; bars[0, a:marks_[i]] = i
            cmap = ListedColormap([COL[g] for g in ORDER]); norm = BoundaryNorm(np.arange(-0.5, len(ORDER) + 0.5), len(ORDER))
        else:
            a = np.where((A.obs["ct"] == "M").values)[0].size
            bars[0, :a] = 0; bars[0, a:] = 1
            cmap = ListedColormap([COL["M"], COL["MyoEpi"]]); norm = BoundaryNorm(np.arange(-0.5, 2.5), 2)
        axb.imshow(bars, aspect="auto", cmap=cmap, norm=norm, extent=(0, nC, 0, 1))
        axb.set_xticks([]); axb.set_yticks([])
        for sp in axb.spines.values(): sp.set_visible(False)
        im = axh.imshow(Z, aspect="auto", cmap="magma", vmin=0, vmax=1, interpolation="nearest")
        axh.set_xticks([]); axh.set_yticks([])
        for sp in axh.spines.values(): sp.set_visible(False)
        axn.set_xlim(0, nC); axn.set_ylim(0, 1); axn.axis("off")
        if axb is axb1:
            for i, g in enumerate(ORDER):
                a = 0 if i == 0 else marks_[i - 1]; b = marks_[i]
                axb.text((a + b) / 2, 0.5, g, ha="center", va="center", fontsize=4.8, weight="bold", color="white")
                if i > 0: axb.axvline(a, color="white", lw=0.7); axh.axvline(a, color="white", lw=0.7)
                axn.text((a + b) / 2, 0.0, "n=%s" % format(total_[i], ","), ha="center", va="bottom", fontsize=3.6)
        else:
            a = np.where((A.obs["ct"] == "M").values)[0].size
            axb.text(a / 2, 0.5, "M", ha="center", va="center", fontsize=4.8, weight="bold", color="white")
            axb.text(a + (nC - a) / 2, 0.5, "MyoEpi", ha="center", va="center", fontsize=4.8, weight="bold", color="white")
            axb.axvline(a, color="white", lw=0.7); axh.axvline(a, color="white", lw=0.7)
            axn.text(a / 2, 0.0, "n=45", ha="center", va="bottom", fontsize=3.6)
            axn.text(a + (nC - a) / 2, 0.0, "n=65", ha="center", va="bottom", fontsize=3.6)
        axb.set_title(label, fontsize=5.6, pad=3)
    axl.set_xlim(0, 1); axl.set_ylim(Z1.shape[0] - 0.5, -0.5); axl.axis("off")
    for i, r in enumerate(P.itertuples()):
        axl.text(1.0, i, "%s (%s)" % (r.sym, r.gid), ha="right", va="center", fontsize=4.3)
    pos = axh2.get_position(); cw, ch = 0.011, pos.height * 0.18
    axc.set_position([pos.x1 + 0.008, pos.y1 - ch - 0.012, cw, ch])
    cb = fig.colorbar(im, cax=axc); cb.set_label("Relative expression", fontsize=4.6)
    cb.ax.tick_params(labelsize=4.0, length=1.0); cb.outline.set_visible(False)
    p = os.path.join(FIG, "Fig2_marker_heatmap_10groups_inset_preview.png")
    fig.savefig(p, dpi=320, bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\fig1_fig2_v2.log", "w", encoding="utf8", buffering=1)
    if "fig1" in sys.argv: fig1_persample()
    if "fig2" in sys.argv: fig2_inset()
    if len(sys.argv) == 1: fig1_persample(); fig2_inset()
