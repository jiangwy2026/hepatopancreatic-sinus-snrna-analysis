# -*- coding: utf-8 -*-
"""Fig. 6b: PH substate composition with corrected substate names."""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

ROOT = Path(r"D:\桌面\single cell")
DATA = ROOT / "analysis_output" / "Fig5_dynamic_subpop_subcluster_composition.csv"
OUTDIR = ROOT / "analysis_output" / "figures_pub" / "2026-10-21_Fig6b"
OUTDIR.mkdir(parents=True, exist_ok=True)

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

d = pd.read_csv(DATA)
d = d[d["population"] == "PH"].sort_values("subcluster")
assert d["subcluster"].tolist() == [0, 1, 2, 3]
assert np.isclose(d["pct_control"].sum(), 100.0)
assert np.isclose(d["pct_hemospasia"].sum(), 100.0)

labels = {
    0: "s0 ECM/ribosome-high",
    1: "s1 progenitor-like",
    2: "s2 immune/ECM-primed",
    3: "s3 TGM1-low",
}
colors = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#54A24B",
}

groups = ["Control", "Hemospasia"]
xpos = np.array([0.0, 0.78])
bottom = np.zeros(2)

fig, ax = plt.subplots(figsize=(4.55, 3.15))
for _, row in d.iterrows():
    sid = int(row["subcluster"])
    vals = np.array([row["pct_control"], row["pct_hemospasia"]])
    ns = np.array([int(row["n_control"]), int(row["n_hemospasia"])])
    ax.bar(
        xpos,
        vals,
        bottom=bottom,
        color=colors[sid],
        width=0.60,
        edgecolor="white",
        linewidth=0.8,
        label=labels[sid],
    )
    for i, (value, n, base) in enumerate(zip(vals, ns, bottom)):
        if value >= 6:
            ax.text(
                xpos[i],
                base + value / 2,
                f"s{sid}  {value:.1f}%\n(n={n})",
                ha="center",
                va="center",
                fontsize=6.3,
                color="white" if sid != 3 else "#1f1f1f",
                linespacing=0.95,
            )
    bottom += vals

ax.set_xticks(xpos)
ax.set_xticklabels(groups)
ax.set_xlim(-0.5, 1.30)
ax.set_ylim(0, 100)
ax.set_yticks(np.arange(0, 101, 20))
ax.set_ylabel("PH cells (%)")
ax.set_title("PH substate composition", pad=7)
ax.text(xpos[0], -8, f"n = {int(d['n_control'].sum())}", ha="center", va="top", fontsize=7)
ax.text(
    xpos[1],
    -8,
    f"n = {int(d['n_hemospasia'].sum())}  (2.30×)",
    ha="center",
    va="top",
    fontsize=7,
)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.tick_params(axis="both", width=0.6, length=3)

handles = [Patch(facecolor=colors[k], edgecolor="none", label=labels[k]) for k in sorted(labels)]
ax.legend(
    handles=handles,
    frameon=False,
    fontsize=6.2,
    loc="center left",
    bbox_to_anchor=(1.02, 0.50),
    handlelength=1.0,
    handletextpad=0.6,
    borderaxespad=0.0,
)

fig.tight_layout()
stem = OUTDIR / "Fig6b_PH_substate_composition"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
fig.savefig(
    stem.with_suffix(".tiff"),
    format="tiff",
    dpi=1200,
    bbox_inches="tight",
    pil_kwargs={"compression": "tiff_lzw"},
)
plt.close(fig)
print(stem)
