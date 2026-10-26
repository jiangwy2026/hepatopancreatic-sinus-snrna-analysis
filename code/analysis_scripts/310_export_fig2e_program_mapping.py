# -*- coding: utf-8 -*-
"""Export Figure 2E: LyO and hemolymph prohemocyte program mapping."""

from pathlib import Path
import collections
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig1_Fig2_final" / "Figure2"
OUT.mkdir(parents=True, exist_ok=True)
A_PATH = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
ELIFE_DIR = ROOT / "analysis_output" / "elife_ref"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.65,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})
ORDER = ["PH", "GH", "TS", "MH2", "MyoEpi", "M", "F", "R", "B", "PM-like"]
PHC = "#B2182B"; GREY = "#8DA0CB"
a = ad.read_h5ad(A_PATH)
T = np.asarray(a.obsm["X_tsne"])
ct = a.obs["cell_type"].astype(str).values

# Combined hemolymph PH1 + PH2 program, using the same mapping logic as the prior analysis.
vn = set(map(str, a.var_names))
sym = collections.defaultdict(list)
for i, s in enumerate(a.var["sym"].astype(str)):
    sym[s.lower()].append(a.var_names[i])

def map_set(fn, n=50):
    d = pd.read_csv(ELIFE_DIR / fn).sort_values("Log2FC", ascending=False).head(n)
    out = []
    for x in d["Gene Name"].astype(str):
        v = x if x in vn else (sym[x.lower()][0] if len(sym.get(x.lower(), [])) == 1 else None)
        if v and v not in out:
            out.append(v)
    return out

p1 = map_set("elife-80127-supp3-v2.csv")
p2 = map_set("elife-80127-supp6-v2.csv")
genes = list(dict.fromkeys(p1 + p2))
idx = [list(a.var_names).index(g) for g in genes]
N = a.n_vars; k = len(idx)
blood = np.zeros(a.n_obs, dtype=np.float32)
Xall = a.X
for st in range(0, a.n_obs, 400):
    en = min(a.n_obs, st + 400)
    X = np.asarray(Xall[st:en].todense(), dtype=np.float32) if hasattr(Xall[st:en], "todense") else np.asarray(Xall[st:en], dtype=np.float32)
    rank = np.argsort(-X, axis=1, kind="stable")
    r = np.empty_like(rank)
    r[np.arange(en-st)[:,None], rank] = np.arange(1, N+1)
    S = r[:, idx].sum(1)
    blood[st:en] = 1 - (S - k*(k+1)/2) / (k*(N-k))

scores = [("LyO progenitor program", a.obs["ucell_LyO_Progenitor"].values),
          ("Hemolymph prohemocyte program", blood)]

fig, axes = plt.subplots(2, 2, figsize=(8.1, 6.7), gridspec_kw={"width_ratios":[1.25,1.0], "hspace":0.34, "wspace":0.34})
for row, (title, score) in enumerate(scores):
    ax = axes[row,0]
    order_idx = np.argsort(score)
    sca = ax.scatter(T[order_idx,0], T[order_idx,1], c=score[order_idx], cmap="viridis", s=2.0, linewidths=0, rasterized=True,
                     vmin=np.quantile(score,0.02), vmax=np.quantile(score,0.98))
    for cell in ORDER:
        m = ct == cell
        x, y = np.median(T[m,0]), np.median(T[m,1])
        ax.text(x, y, cell, fontsize=6.2, ha="center", va="center", color="#1e293b",
                bbox=dict(boxstyle="round,pad=0.13", fc="white", ec="none", alpha=0.78), zorder=4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("tSNE"); ax.set_ylabel("tSNE")
    ax.set_title(title, pad=4)
    ax.set_box_aspect(1)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    cb = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("UCell score", fontsize=7); cb.ax.tick_params(labelsize=6.2)

    ax2 = axes[row,1]
    means = [score[ct == cell].mean() for cell in ORDER]
    pct = [np.mean(score[ct == cell] > np.median(score))*100 for cell in ORDER]
    ax2.hlines(range(len(ORDER)), min(means)*0.99, means, color="#DDDDDD", lw=1.0)
    ax2.scatter(means, range(len(ORDER)), s=[20+p*1.7 for p in pct],
                c=[PHC if cell=="PH" else GREY for cell in ORDER], edgecolors="#333", linewidths=0.4, zorder=3)
    for i,(m,p) in enumerate(zip(means,pct)):
        ax2.text(m+0.002, i, f"{m:.3f} ({p:.0f}%)", va="center", fontsize=5.9)
    ax2.set_yticks(range(len(ORDER))); ax2.set_yticklabels(ORDER); ax2.invert_yaxis()
    ax2.set_xlim(min(means)*0.985, max(means)*1.06)
    ax2.set_xlabel("Mean UCell score")
    ax2.set_title(title+" summary", pad=4)
    for sp in ("top","right"): ax2.spines[sp].set_visible(False)

fig.tight_layout()
stem = OUT / "Fig2E_program_mapping"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
fig.savefig(stem.with_suffix(".tiff"), format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)
print(stem)
