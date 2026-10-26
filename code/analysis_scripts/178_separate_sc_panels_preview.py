# -*- coding: utf-8 -*-
"""178_separate_sc_panels_preview.py - render each new-name panel as a standalone preview."""
import os, sys, importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\桌面\single cell\analysis_output"
OUT = os.path.join(BASE, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(OUT, exist_ok=True)

spec = importlib.util.spec_from_file_location("sc176", os.path.join(BASE, "scripts", "176_sc_summary_new_names_preview.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def save(fig, name, dpi=400):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white", pad_inches=0.04)
    print("saved", p)
    plt.close(fig)

def hide_letters(ax, letters=("a","b","c","d")):
    for t in ax.texts:
        if t.get_text() in letters:
            t.set_visible(False)

# 1 composition
fig = plt.figure(figsize=(4.8, 3.35), facecolor="white")
ax = fig.add_axes([0.34, 0.17, 0.62, 0.72])
m.draw_composition(ax)
save(fig, "01_composition_preview.png")

# 2 pseudotime
fig = plt.figure(figsize=(4.1, 4.1), facecolor="white")
ax = fig.add_axes([0.12, 0.12, 0.72, 0.76])
m.fig = fig
m.draw_pseudotime(ax)
hide_letters(ax)
save(fig, "02_pseudotime_preview.png")

# 3 velocity
fig = plt.figure(figsize=(4.4, 4.8), facecolor="white")
gs = fig.add_gridspec(2, 1, height_ratios=[2.0, 1.0], hspace=0.65,
                      left=0.15, right=0.97, top=0.91, bottom=0.10)
ax1 = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
m.fig = fig
m.draw_velocity(ax1, ax2)
hide_letters(ax1); hide_letters(ax2)
save(fig, "03_velocity_summary_preview.png")

# 4 mechanism
fig = plt.figure(figsize=(4.8, 5.1), facecolor="white")
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.55,
                      left=0.13, right=0.97, top=0.93, bottom=0.09)
ax1 = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
m.draw_mechanism(ax1, ax2)
save(fig, "04_ccnd2_TGM_mechanism_preview.png")
