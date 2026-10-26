# -*- coding: utf-8 -*-
"""219_dynamic_subcluster_pseudotime.py
Detailed pseudotime view of eligible dynamic substates using the existing
global velocity pseudotime. Generates local UMAP and ridge-density previews.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, mannwhitneyu

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
SUB_DIR = os.path.join(OUT, "dynamic_subpop_subclusters")
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
sc.settings.verbosity = 0

GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
PARENTS = ["PH", "MH2", "R", "PM-like"]
PALETTE = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#9C755F"]
SHORT = {
    "PH_s0": "s0 ECM-high", "PH_s1": "s1 progenitor", "PH_s2": "s2 immune-primed", "PH_s3": "s3 TGM1+",
    "MH2_s0": "s0 SDC2/TGM1+", "MH2_s1": "s1 ECM-remodeling", "MH2_s2": "s2 SSPO+", "MH2_s3": "s3 activated macrophage",
    "R_s0": "s0 transitional", "R_s2": "s2 neprilysin+", "R_s3": "s3 SLC5A1/NAAT1+",
    "PM-like_s0": "s0 core", "PM-like_s2": "s2 stress-high core",
}


def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0: return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order); prev = 1.0
    for rank, i in enumerate(order[::-1], 1):
        val = min(prev, p[i] * m / (m - rank + 1))
        out[i] = val; prev = val
    return out


def read_subclusters():
    frames = []
    for f in glob.glob(os.path.join(SUB_DIR, "*_subcluster_umap.csv")):
        grp = os.path.basename(f).split("_")[0]
        d = pd.read_csv(f, index_col=0)
        d["parent"] = grp
        d["sub"] = grp + "_s" + d["subcluster"].astype(str)
        frames.append(d)
    return pd.concat(frames)


def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"), backed="r")
    obs = A.obs.copy(); A.file.close()
    obs["parent"] = obs["vendor_cluster"].astype(str).map(GROUPS)
    obs["group2"] = obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    S = read_subclusters()
    obs = obs.join(S[["sub"]], how="left")
    T = pd.read_csv(os.path.join(OUT, "Fig5_dynamic_subcluster_annotation_table.csv"))
    E = T[T.analysis_status == "eligible"].copy()

    rows = []
    for r in E.itertuples():
        d = obs[obs["sub"].astype(str) == r.key]
        for g in ("Control", "Hemospasia"):
            x = d.loc[d.group2 == g, "velocity_pseudotime"].dropna().astype(float)
            rows.append({"parent": r.parent, "subcluster": r.subcluster, "key": r.key,
                         "group": g, "n": int(len(x)), "median_pseudotime": float(x.median()) if len(x) else np.nan,
                         "q25": float(x.quantile(.25)) if len(x) else np.nan,
                         "q75": float(x.quantile(.75)) if len(x) else np.nan})
    R = pd.DataFrame(rows)
    pvals = []
    keys = E.key.tolist()
    for key in keys:
        c = obs[(obs["sub"].astype(str) == key) & (obs.group2 == "Control")]["velocity_pseudotime"].dropna().astype(float)
        t = obs[(obs["sub"].astype(str) == key) & (obs.group2 == "Hemospasia")]["velocity_pseudotime"].dropna().astype(float)
        if len(c) >= 5 and len(t) >= 5:
            pvals.append(mannwhitneyu(c, t, alternative="two-sided").pvalue)
        else:
            pvals.append(np.nan)
    qvals = bh(pvals)
    P = pd.DataFrame({"key": keys, "mannwhitney_p": pvals, "fdr": qvals})
    R = R.merge(P, on="key", how="left")
    E = E.merge(P, on="key", how="left")
    sig_keys = set(E.loc[E.fdr < 0.05, "key"])
    R.to_csv(os.path.join(OUT, "Fig5_dynamic_subcluster_pseudotime_summary.csv"),
             index=False, encoding="utf-8-sig")

    # Ridge-density preview
    fig = plt.figure(figsize=(7.2, 5.0), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.98, top=0.88, bottom=0.12,
                          wspace=0.30, hspace=0.48)
    for pi, parent in enumerate(PARENTS):
        ax = fig.add_subplot(gs[pi // 2, pi % 2])
        ek = E[E.parent == parent].key.tolist()
        dpar = obs[(obs.parent == parent) & obs["sub"].astype(str).isin(ek)].copy()
        xmin, xmax = dpar.velocity_pseudotime.quantile([0.01, 0.99]).astype(float).values
        pad = 0.08 * (xmax - xmin if xmax > xmin else 1)
        xx = np.linspace(xmin - pad, xmax + pad, 220)
        for row, key in enumerate(ek[::-1]):
            col = PALETTE[row % len(PALETTE)]
            dsub = dpar[dpar["sub"].astype(str) == key]
            for g, ls, alpha in [("Control", (0, (3, 2)), 0.75), ("Hemospasia", "-", 1.0)]:
                x = dsub.loc[dsub.group2 == g, "velocity_pseudotime"].dropna().astype(float).values
                if len(x) >= 5 and x.std() > 1e-9:
                    y = gaussian_kde(x, bw_method="scott")(xx)
                    y = y / max(y.max(), 1e-9) * 0.80
                    ax.plot(xx, row + y, color=col, lw=0.9, ls=ls, alpha=alpha)
                    ax.fill_between(xx, row, row + y, color=col, alpha=0.08 if g == "Hemospasia" else 0)
                    ax.plot([np.median(x), np.median(x)], [row - 0.05, row + 0.10], color=col,
                            lw=1.1 if g == "Hemospasia" else 0.7, ls=ls)
        ax.set_yticks(np.arange(len(ek)))
        ax.set_yticklabels([SHORT.get(k, k) + (" *" if k in sig_keys else "") for k in ek[::-1]], fontsize=4.9)
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(-0.35, len(ek) - 0.35)
        ax.set_xlabel("Velocity pseudotime", fontsize=5.2)
        ax.tick_params(axis="x", labelsize=4.9, length=1.3, width=0.4)
        ax.tick_params(axis="y", length=0, pad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="x", color="#EFEFEF", lw=0.3, zorder=0)
        ax.set_title(f"{parent}", fontsize=6.2, fontweight="bold", pad=3)
    fig.text(0.50, 0.96, "Pseudotime distribution of eligible dynamic substates",
             ha="center", va="top", fontsize=7.2)
    fig.text(0.50, 0.025, "Solid = Hemospasia; dashed = Control; vertical ticks = medians. "
             "* FDR < 0.05 for group difference (cell-level Mann–Whitney). "
             "Global velocity pseudotime is a transcriptional ordering, not lineage tracing.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out1 = os.path.join(PARTS, "Fig5_dynamic_subcluster_pseudotime_density_preview.png")
    fig.savefig(out1, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out1)
    plt.close(fig)

    # Local UMAP coloured by global pseudotime
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), facecolor="white")
    for ax, parent in zip(axes.flat, PARENTS):
        d = S[S.parent == parent].copy()
        d = d.reindex(obs.index.intersection(d.index))
        pt = obs.loc[d.index, "velocity_pseudotime"].astype(float)
        d["pseudotime"] = pt.rank(pct=True).values
        sca = ax.scatter(d["umap1"], d["umap2"], c=d["pseudotime"], s=2.2, cmap="viridis",
                         linewidths=0, rasterized=True, alpha=0.9, vmin=0, vmax=1)
        for key in E[E.parent == parent].key.tolist():
            dsub = d[d["sub"] == key]
            if len(dsub) >= 3:
                ax.text(dsub.umap1.median(), dsub.umap2.median(), key.split("_")[1],
                        fontsize=4.8, ha="center", va="center", color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.10", fc="#333333", ec="none", alpha=0.65))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.35); sp.set_color("#777777")
        ax.set_title(parent, fontsize=6.2, fontweight="bold", pad=3)
    fig.suptitle("Global velocity pseudotime projected onto substate UMAPs", y=0.98, fontsize=7.2)
    cax = fig.add_axes([0.92, 0.20, 0.015, 0.55])
    cb = fig.colorbar(sca, cax=cax); cb.set_label("Within-population pseudotime percentile", fontsize=5.2, labelpad=3)
    cb.ax.tick_params(labelsize=4.9, length=1.3, width=0.4); cb.outline.set_linewidth(0.4)
    fig.text(0.5, 0.02, "Only eligible substates are labelled; small excluded substates remain visible as cells.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out2 = os.path.join(PARTS, "Fig5_dynamic_subcluster_pseudotime_umap_preview.png")
    fig.savefig(out2, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out2)
    plt.close(fig)


if __name__ == "__main__":
    main()
