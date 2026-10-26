# -*- coding: utf-8 -*-
"""179_composition_stacked_preview.py - standalone stacked composition preview."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import numpy as np

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)

OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
COL = {"PH": "#B2182B", "GH": "#E08214", "TS": "#D6604D", "MH2": "#8C510A",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "M": "#2166AC",
       "MyoEpi": "#92C5DE", "PM-like": "#7B3294"}
SHORT = {"PM-like": "PM-like", "MyoEpi": "MyoEpi"}

def main():
    comp = pd.read_csv(os.path.join(OUT, "composition_new_names_10groups.csv"))
    chi = json.load(open(os.path.join(OUT, "composition_new_names_chi2.json")))
    order = comp["group"].tolist()
    fig = plt.figure(figsize=(4.8, 3.55), facecolor="white")
    ax = fig.add_axes([0.14, 0.37, 0.82, 0.50])
    yvals = [1.0, 0.0]
    ylabels = ["Control (0 h)\nn = %s" % format(chi["n_blank"], ","),
               "Hemospasia (6 h)\nn = %s" % format(chi["n_treated"], ",")]
    left = np.zeros(2)
    for g in order:
        vals = np.array([comp.loc[comp.group == g, "pct_blank"].iloc[0],
                         comp.loc[comp.group == g, "pct_treated"].iloc[0]], dtype=float)
        ax.barh(yvals, vals, left=left, color=COL[g], height=0.52,
                edgecolor="white", linewidth=0.55, zorder=2)
        for y, l, v in zip(yvals, left, vals):
            if v >= 4.0:
                ax.text(l + v / 2, y, "%.1f" % v, ha="center", va="center",
                        fontsize=5.0, color="white" if g not in ("B", "MyoEpi") else "#222222",
                        zorder=3)
        left += vals
    ax.set_yticks(yvals); ax.set_yticklabels(ylabels, fontsize=5.8)
    ax.set_xlim(0, 100); ax.set_ylim(-0.45, 1.45)
    ax.set_xlabel("Cells (%)", fontsize=6.0)
    ax.tick_params(axis="x", labelsize=5.2, length=1.7, width=0.5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E4E4E4", lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Cell composition", fontsize=7.3, pad=6)

    # Matplotlib fills legend columns first; interleave for row-wise reading.
    half = (len(order) + 1) // 2
    legend_order = []
    for i in range(half):
        legend_order.append(order[i])
        if i + half < len(order):
            legend_order.append(order[i + half])
    handles = [Patch(facecolor=COL[g], edgecolor="white", label="%s (%s)" % (g, format(int(comp.loc[comp.group == g, "n_total"].iloc[0]), ",")))
               for g in legend_order]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.205),
               ncol=5, frameon=False, fontsize=5.1, handlelength=0.9,
               handletextpad=0.28, columnspacing=0.85, labelspacing=0.35)

    p = os.path.join(PREV, "01b_composition_stacked_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)

if __name__ == "__main__":
    main()
