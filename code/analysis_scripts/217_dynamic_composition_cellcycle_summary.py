# -*- coding: utf-8 -*-
"""217_dynamic_composition_cellcycle_summary.py
Composition effect sizes for all ten populations and cell-cycle phase
composition for the five dynamic populations, split by group."""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-20", "第二张图_抽血处理细胞变化")
os.makedirs(PARTS, exist_ok=True); os.makedirs(DELIV, exist_ok=True)
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
DYN = ["PH", "MH2", "R", "MyoEpi", "PM-like"]
UP = "#B2182B"; DOWN = "#2166AC"; INK = "#222222"
PHCOL = {"G1": "#BFBFBF", "S": "#4393C3", "G2M": "#B2182B"}
S_GENES = ["PCNA", "MCM6", "CDC6", "ORC1", "ORC4", "ORC5", "POLA2", "PRIM1",
           "FEN1", "TYMS", "DHFR", "HIRA", "CHAF1A"]
G2M_GENES = ["CCNB3", "AURKA", "TOP2A", "TPX2", "KIF20A", "KIF23", "RACGAP1",
             "NDC80", "INCENP", "PRC1", "ASPM", "Pka-C1", "Patronin"]


def main():
    comp = pd.read_csv(os.path.join(OUT, "composition_new_names_10groups.csv"))
    comp["log2_fold"] = np.log2(comp["fold_treated_vs_blank"])
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.obs["ct"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    A.obs["group2"] = A.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    sym_to_gid = {}
    for gid, sym in zip(A.var_names.astype(str), A.var["sym"].astype(str)):
        sym_to_gid.setdefault(sym, gid)
    s_ids = [sym_to_gid[g] for g in S_GENES if g in sym_to_gid]
    g2m_ids = [sym_to_gid[g] for g in G2M_GENES if g in sym_to_gid]
    print("cell-cycle genes:", len(s_ids), "S /", len(g2m_ids), "G2M")
    sc.tl.score_genes_cell_cycle(A, s_genes=s_ids, g2m_genes=g2m_ids, random_state=0)
    cc_rows = []
    for grp in DYN:
        sub = A[A.obs["ct"].astype(str) == grp]
        tab = pd.crosstab(sub.obs["group2"].astype(str), sub.obs["phase"].astype(str))
        tab = tab.reindex(columns=["G1", "S", "G2M"], fill_value=0)
        for group, r in tab.iterrows():
            total = int(r.sum())
            row = {"population": grp, "group": group, "n": total}
            for ph in ["G1", "S", "G2M"]:
                row[f"pct_{ph}"] = 100 * float(r[ph]) / total if total else np.nan
            cc_rows.append(row)
    CC = pd.DataFrame(cc_rows)
    CC.to_csv(os.path.join(OUT, "Fig5_dynamic_cellcycle_by_group.csv"), index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 4.2), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.42], left=0.075, right=0.975,
                          top=0.87, bottom=0.17, wspace=0.42)
    # composition effect sizes
    ax1 = fig.add_subplot(gs[0, 0])
    d = comp.sort_values("log2_fold")
    y = np.arange(len(d))
    colors = [UP if v > 0 else DOWN for v in d.log2_fold]
    ax1.barh(y, d.log2_fold, color=colors, height=0.62, edgecolor="white", linewidth=0.4)
    ax1.axvline(0, color="#333333", lw=0.7)
    ax1.set_yticks(y)
    ax1.set_yticklabels([f"{r.group} ({int(r.n_total):,})" for r in d.itertuples()], fontsize=5.2)
    ax1.set_xlabel("log$_2$ fold change (Hemospasia / Control)", fontsize=5.5)
    ax1.tick_params(axis="x", labelsize=5.0, length=1.5, width=0.4)
    ax1.tick_params(axis="y", length=0, pad=2)
    for i, (v, fold) in enumerate(zip(d.log2_fold, d.fold_treated_vs_blank)):
        ax1.text(v + (0.04 if v >= 0 else -0.04), i, f"{fold:.2f}×",
                 ha="left" if v >= 0 else "right", va="center", fontsize=4.9, color=INK)
    ax1.set_xlim(-2.5, 1.7)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.grid(axis="x", color="#EEEEEE", lw=0.35, zorder=0)
    ax1.set_title("Cell composition change", fontsize=6.2, pad=4)

    # cell-cycle composition
    ax2 = fig.add_subplot(gs[0, 1])
    order = []
    for grp in DYN:
        order += [(grp, "Control"), (grp, "Hemospasia")]
    ypos = []
    cur = 0
    for grp in DYN:
        ypos += [cur + 0.18, cur - 0.18]
        cur -= 1.0
    for (grp, group), yy in zip(order, ypos):
        r = CC[(CC.population == grp) & (CC.group == group)].iloc[0]
        left = 0.0
        for ph in ["G1", "S", "G2M"]:
            v = float(r[f"pct_{ph}"])
            ax2.barh(yy, v, left=left, color=PHCOL[ph], height=0.30,
                     edgecolor="white", linewidth=0.35)
            if v >= 10:
                ax2.text(left + v / 2, yy, f"{v:.0f}", ha="center", va="center",
                         fontsize=4.6, color="#333333" if ph == "G1" else "white")
            left += v
    ax2.set_yticks(ypos)
    labels = []
    for grp, group in order:
        n = int(CC[(CC.population == grp) & (CC.group == group)].iloc[0].n)
        labels.append(f"{grp} {group} (n={n})")
    ax2.set_yticklabels(labels, fontsize=4.7)
    ax2.set_xlim(0, 100); ax2.set_xlabel("Cells per phase (%)", fontsize=5.5)
    ax2.tick_params(axis="x", labelsize=5.0, length=1.5, width=0.4)
    ax2.tick_params(axis="y", length=0, pad=2)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.grid(axis="x", color="#EEEEEE", lw=0.35, zorder=0)
    ax2.legend(handles=[Patch(color=PHCOL[p], label=p) for p in ["G1", "S", "G2M"]],
               loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=3, frameon=False,
               fontsize=5.0, handlelength=0.8, columnspacing=1.0)
    ax2.set_title("Cell-cycle phase by group", fontsize=6.2, pad=4)

    fig.text(0.50, 0.97, "Composition change and cell-cycle state", ha="center", va="top",
             fontsize=7.2)
    fig.text(0.50, 0.035, "Cell-cycle scores use 13 S and 13 G2M genes. "
             "Cell-level composition is descriptive; no biological replicates.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out = os.path.join(PARTS, "Fig5_dynamic_composition_cellcycle_preview.png")
    fig.savefig(out, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(os.path.join(DELIV, "Fig5_composition_cellcycle.svg"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(os.path.join(DELIV, "Fig5_composition_cellcycle.tiff"), dpi=1200, bbox_inches="tight", pad_inches=0.04, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    print("saved", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
