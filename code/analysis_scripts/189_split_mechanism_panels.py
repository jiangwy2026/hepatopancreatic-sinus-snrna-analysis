# -*- coding: utf-8 -*-
"""189_split_mechanism_panels.py - standalone volcano and violin previews."""
import os, sys, importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(BASE, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
spec = importlib.util.spec_from_file_location("m188", os.path.join(BASE, "scripts", "188_ccnd2_TGM_mechanism_preview.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def save(fig, name):
    p = os.path.join(PREV, name)
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)

# volcano only
fig = plt.figure(figsize=(4.35, 4.05), facecolor="white")
ax = fig.add_axes([0.13, 0.14, 0.82, 0.79])
m.draw_volcano(ax)
save(fig, "04a_volcano_earlyPH_preview.png")

# violin only
fig = plt.figure(figsize=(4.35, 4.05), facecolor="white")
ax = fig.add_axes([0.14, 0.14, 0.82, 0.79])
m.draw_violin(ax)
ax.set_title("Early PH: ccnd2 and hemocyte TGM", fontsize=6.8, pad=4)
save(fig, "04b_violin_ccnd2_TGM_preview.png")
