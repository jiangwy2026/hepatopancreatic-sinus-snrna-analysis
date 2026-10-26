# -*- coding: utf-8 -*-
"""Supp13 PH substate bubble plot with colorbar and size legend aligned in one column."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Supp13_bubble_same_column"
OUT.mkdir(parents=True, exist_ok=True)
H5 = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
PH = ROOT / "analysis_output" / "dynamic_subpop_subclusters" / "PH_subcluster_umap.csv"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6,
    "savefig.facecolor": "white",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

a = ad.read_h5ad(H5, backed="r")
ph = pd.read_csv(PH)
ph.columns = ["cell", "group", "sub", "u1", "u2"]
sub = ph.set_index("cell")["sub"]
keep = np.isin(a.obs_names, sub.index)
A = a[keep]
s = sub.reindex(A.obs_names).values
genes = [
    ("s0", "PEN-3", "ncbi_113808997"), ("s0", "LOC113819543", "ncbi_113819543"),
    ("s0", "LOC113819837", "ncbi_113819837"), ("s0", "LOC113830563", "ncbi_113830563"),
    ("s0", "LOC113802679", "ncbi_113802679"), ("s1", "LOC113808957", "ncbi_113808957"),
    ("s1", "LOC113808628", "ncbi_113808628"), ("s1", "LOC113800422", "ncbi_113800422"),
    ("s1", "Igfbp7", "ncbi_113811702"), ("s1", "ccnd2", "ncbi_113814652"),
    ("s2", "LOC113810024", "ncbi_113810024"), ("s2", "THBS4", "ncbi_113823118"),
    ("s2", "TGM1", "ncbi_113823945"), ("s3", "LOC113818748", "gene-LOC113818748"),
    ("s3", "LOC113823090", "gene-LOC113823090"), ("s3", "LOC113823089", "gene-LOC113823089"),
    ("context", "PPAF2", "ncbi_113812814"), ("context", "Grik2", "ncbi_113823802"),
]
rows = []
for grp, name, key in genes:
    hit = key if key in A.var_names else key.replace("gene-", "")
    if hit not in A.var_names:
        print("missing", name, key)
        continue
    x = np.asarray(A[:, hit].layers["lognorm"].todense()).ravel()
    for s_ in [0, 1, 2, 3]:
        m = s == s_
        rows.append({"group": grp, "gene": name, "sub": s_, "mean": float(x[m].mean()), "pct": float((x[m] > 0).mean() * 100)})
d = pd.DataFrame(rows)
d.to_csv(OUT / "Supp13_PH_substate_markers_v3.csv", index=False, encoding="utf-8-sig")
order = d.drop_duplicates("gene").gene.tolist()

X, Y, SZ, CC = [], [], [], []
for i, g in enumerate(order):
    for j in [0, 1, 2, 3]:
        r = d[(d.gene == g) & (d["sub"] == j)]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        X.append(j); Y.append(i); SZ.append(max(r.pct * 1.1, 2)); CC.append(r["mean"])

fig = plt.figure(figsize=(5.9, 5.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.25], height_ratios=[1.0, 1.0],
                      left=0.18, right=0.93, bottom=0.11, top=0.92, wspace=0.08, hspace=0.18)
ax = fig.add_subplot(gs[:, 0])
sc = ax.scatter(X, Y, s=SZ, c=CC, cmap="Reds", vmin=0, vmax=d["mean"].max(),
                edgecolors="#666", linewidths=0.25, rasterized=True)
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(["s0\nECM/ribosome", "s1\nPH progenitor", "s2\nImmune/ECM", "s3\nTGM1-low"], fontsize=6.0)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=6.2)
ax.set_xlim(-0.6, 3.6); ax.set_ylim(len(order) - 0.5, -0.5)
ax.set_title("PH substate markers", pad=4)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# Colorbar in the upper half of the shared right column.
cax = fig.add_subplot(gs[0, 1])
cb = fig.colorbar(sc, cax=cax, orientation="vertical")
cb.set_label("Mean expression\n(log1p CP10K)", fontsize=5.8, labelpad=4)
cb.ax.tick_params(labelsize=5.4, length=2)

# Bubble-size legend in the lower half of the same right column.
leg = fig.add_subplot(gs[1, 1])
leg.axis("off")
leg.text(0.0, 0.95, "% expressing", fontsize=6.2, fontweight="bold", ha="left", va="top")
for k, v in enumerate([0, 25, 50, 75]):
    yy = 0.74 - k * 0.19
    leg.scatter(0.18, yy, s=max(v * 1.1, 14), c="#BBBBBB", edgecolors="#666", linewidths=0.4)
    leg.text(0.42, yy, f"{v}%", va="center", fontsize=6.0)
leg.set_xlim(0, 1); leg.set_ylim(0, 1)

stem = OUT / "Supp13_PH_substate_markers_legends_same_column"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
fig.savefig(stem.with_suffix(".tiff"), format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(stem)
