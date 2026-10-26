# -*- coding: utf-8 -*-
"""116_export_cellcycle.py - export the approved cell-cycle bar chart in all four formats."""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ROWS = ["PH", "GH", "MH1", "MH2"]
PHCOL = {"G1": "#BFBFBF", "S": "#4393C3", "G2M": "#B2182B"}

def main():
    comp = pd.read_csv(os.path.join(OUT, "cellcycle_phase_composition.csv"), index_col=0).reindex(ROWS)
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    left = np.zeros(len(ROWS))
    for ph in ("G1", "S", "G2M"):
        v = comp[ph].values
        ax.barh(np.arange(len(ROWS)), v, left=left, color=PHCOL[ph], edgecolor="white", height=0.62, label=ph)
        for k, (l, x) in enumerate(zip(left, v)):
            if x >= 5:
                ax.text(l + x / 2, k, "%.0f" % x, fontsize=6, ha="center", va="center",
                        color="#333333" if ph == "G1" else "white")
        left += v
    ax.set_xlim(0, 100); ax.set_ylim(len(ROWS) - 0.45, -0.55)
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS, fontsize=8, fontweight="bold")
    ax.set_xlabel("Cells per phase (%)", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=6.5, length=1.8, width=0.5); ax.tick_params(axis="y", length=0, pad=3)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False, fontsize=6.5,
              handlelength=0.9, handletextpad=0.3, columnspacing=0.9)
    bp = os.path.join(FIG, "Fig_cellcycle_bar")
    for e, kw in (("_preview.png", dict(dpi=400)), (".pdf", {}), (".svg", {}),
                  (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + e, bbox_inches="tight", pad_inches=0.03, **kw)
        print("  %-40s %6.2f MB" % (os.path.basename(bp + e), os.path.getsize(bp + e) / 1e6))
    plt.close(fig)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\export_cellcycle.log", "w", encoding="utf8", buffering=1)
    main()
