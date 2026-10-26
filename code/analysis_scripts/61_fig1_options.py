import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")
base = {"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"}
A.obs["c10"] = A.obs.vendor_cluster.astype(str).map(base)
A.obs["c8"]  = A.obs["c10"].replace({"NE":"NEU","Neu":"NEU"})
A.obs["c7"]  = A.obs["c8"].replace({"F":"HP","R":"HP","B":"HP"})
FULL = {"PH":"Prohemocyte","GH":"Granular hemocyte","MH1":"Monocytic hemocyte 1","MH2":"Macrophage-like phagocyte",
        "Epi":"Epidermal / cuticle cell","F":"Fibrillar cell","R":"Resorptive cell","B":"Blister-like cell",
        "NE":"Neuroendocrine cell","Neu":"Neural cell","NEU":"Neural / neuroendocrine","HP":"Hepatopancreas cell"}
C10 = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294","F":"#006D2C","R":"#41AB5D","B":"#78C679","NE":"#4393C3","Neu":"#2166AC"}
C8  = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294","F":"#006D2C","R":"#41AB5D","B":"#78C679","NEU":"#2166AC"}
C7  = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294","HP":"#41AB5D","NEU":"#2166AC"}
um = A.obsm["X_umap"]
def panel(ax, key, cols, title, labpos=None):
    for ab in cols:
        m = (A.obs[key] == ab).values
        ax.scatter(um[m,0], um[m,1], s=1.8, c=cols[ab], linewidths=0, rasterized=True,
                   label="%s  %s (%d)" % (ab, FULL[ab], m.sum()))
    for ab in cols:
        m = (A.obs[key] == ab).values
        if m.sum() == 0: continue
        pos = (labpos or {}).get(ab)
        cx, cy = pos if pos else (np.median(um[m,0]), np.median(um[m,1]))
        ax.text(cx, cy, ab, fontsize=6.5, weight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.82))
    ax.set_aspect("equal"); ax.set_title(title, fontsize=8, pad=5)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.tick_params(axis="both", labelsize=5.5, length=2, width=0.6)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=5.4, ncol=1,
              handletextpad=0.3, labelspacing=0.32, markerscale=3.5, frameon=False)
fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.8))
panel(axes[0], "c10", C10, "Option A: 10 populations (current)")
panel(axes[1], "c8",  C8,  "Option B: 8 populations (merge NE + Neu)")
panel(axes[2], "c7",  C7,  "Option C: 7 populations (also merge F + R + B)",
      labpos={"HP": (10.5, 5.6), "B": (5.9, 2.4), "R": (12.3, 2.4)})
plt.tight_layout()
p = os.path.join(figdir, "Fig1_UMAP_options_10_8_7.png")
fig.savefig(p, dpi=200, bbox_inches="tight", pad_inches=0.05); plt.close()
print("ok:", p)
