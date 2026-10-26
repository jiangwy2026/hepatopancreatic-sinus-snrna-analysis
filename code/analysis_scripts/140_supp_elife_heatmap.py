# -*- coding: utf-8 -*-
"""140_supp_elife_heatmap.py - 补充图：eLife 血细胞 6 亚型 marker 集 × 本数据 10 群的标签转移 z-score"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
CL = ["0","1","2","3","4","5","6","7","8","9"]
NAME = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NEU","9":"NEU"}
EN = {"PH1(原血细胞1)": "PH1 (prohemocyte 1)", "PH2(原血细胞2)": "PH2 (prohemocyte 2)",
      "GH1(颗粒1)": "GH1 (granular 1)", "GH2(颗粒2)": "GH2 (granular 2)",
      "MH1(单核1)": "MH1 (monocytic 1)", "MH2(单核2·巨噬样)": "MH2 (macrophage-like)"}

def main():
    E = pd.read_csv(os.path.join(OUT, "elife_label_transfer_z.csv"), index_col=0)
    E.columns = [str(c) for c in E.columns]
    E.index = [EN.get(str(i), str(i)) for i in E.index]
    M = E[CL].astype(float).values
    fig, ax = plt.subplots(figsize=(4.6, 2.3))
    vmax = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax), aspect="auto")
    ax.set_xticks(range(len(CL)))
    ax.set_xticklabels(["%s\n(C%s)" % (NAME[c], c) for c in CL], fontsize=6)
    ax.set_yticks(range(len(E))); ax.set_yticklabels(E.index, fontsize=6.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if abs(M[i, j]) > 1.5:
                ax.text(j, i, "%.1f" % M[i, j], fontsize=5, ha="center", va="center",
                        color="white" if abs(M[i, j]) > 2 else "#333333")
    ax.set_xlabel("Our clusters", fontsize=7)
    ax.set_title("Label transfer: published hemocyte subtypes\n(Yang et al. 2022, eLife; marker match 91.7-98.3%)", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("marker-set enrichment (z)", fontsize=6); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
    for s in ax.spines.values(): s.set_visible(False)
    p = os.path.join(FIG, "FigS_annotation_elife_heatmap_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # 同时导出证据表
    mk = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected_fixed.csv"))
    rows = []
    for c in CL:
        rows.append({"cluster": "C" + c, "our_name": NAME[c], "top_markers": ", ".join(mk[mk.grp == NAME[c]]["gene"].astype(str).tolist()[:5]),
                     "best_eLife_subtype": E[c].idxmax(), "z": round(float(E[c].max()), 2)})
    T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "supp_elife_annotation_table.csv"), index=False)
    print(T.to_string(index=False))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\supp_elife.log", "w", encoding="utf8", buffering=1)
    main()
