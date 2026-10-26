# -*- coding: utf-8 -*-
"""218_dynamic_subcluster_annotation_summary.py
Assign cautious functional labels to eligible subclusters, quantify their
composition changes, and export annotation/change summary tables + previews.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
SUB_DIR = os.path.join(OUT, "dynamic_subpop_subclusters")
os.makedirs(PARTS, exist_ok=True)
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
ORDER_PARENT = ["PH", "MH2", "R", "PM-like", "MyoEpi"]
PROGRAMS = {
    "PH/progenitor": ["PPAF2", "SSPO", "Grik2", "ccnd2", "TGM1", "Syne1", "slo", "SERPINB1"],
    "Macrophage/immune": ["TIE1", "ITGA2B", "PTPRA", "SDC2", "Tl", "AGAP007347"],
    "Absorptive/metabolic": ["ACE", "Oacyl", "HSPG2", "mesh", "hdc", "NAAT1", "SLC5A1", "Npc1"],
    "ECM/remodeling": ["THBS4", "HSPG2", "Tinagl1", "Angpt4", "SDC2"],
    "Stress": ["hsc71", "Hspa14", "HSP90", "HSPA", "DnaJ"],
    "PM-like/cuticle": ["Fas2", "LOC113809538", "LOC113809539", "LOC113809540", "GPx"],
    "G2M": ["CCNB3", "AURKA", "TOP2A", "TPX2", "KIF20A", "KIF23", "RACGAP1", "NDC80", "INCENP", "PRC1", "ASPM", "Pka-C1", "Patronin"],
}
LABELS = {
    "PH_s0": "ECM/ribosome-high PH state",
    "PH_s1": "PH progenitor-like state",
    "PH_s2": "Immune-primed PH state",
    "PH_s3": "TGM1+ progenitor-like state",
    "MH2_s0": "SDC2/TGM1+ immune-ECM state",
    "MH2_s1": "THBS4/Tinagl1+ ECM-remodeling state",
    "MH2_s2": "SSPO+ immune-ECM state",
    "MH2_s3": "Angpt4/SDC2+ activated macrophage-like state",
    "R_s0": "PEN-3+ transitional R state",
    "R_s1": "PEN-3+ R state (too few cells)",
    "R_s2": "Neprilysin+ absorptive R state",
    "R_s3": "SLC5A1/NAAT1+ absorptive R state",
    "R_s4": "R state (too few cells)",
    "PM-like_s0": "PM-like core/cuticle state",
    "PM-like_s2": "PM-like stress-high core state",
}
LABEL_ORDER = ["PH_s0", "PH_s1", "PH_s2", "PH_s3",
               "MH2_s0", "MH2_s1", "MH2_s2", "MH2_s3",
               "R_s0", "R_s2", "R_s3",
               "PM-like_s0", "PM-like_s2"]


def main():
    comp = pd.read_csv(os.path.join(OUT, "Fig5_dynamic_subpop_subcluster_composition.csv"))
    comp["key"] = comp["population"].astype(str) + "_s" + comp["subcluster"].astype(str)
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.obs["parent"] = A.obs["vendor_cluster"].astype(str).map(GROUPS)
    A.obs["grp2"] = A.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    subcol = pd.Series(index=A.obs_names, dtype=object)
    for f in glob.glob(os.path.join(SUB_DIR, "*_subcluster_umap.csv")):
        grp = os.path.basename(f).split("_")[0]
        d = pd.read_csv(f, index_col=0)
        subcol.loc[d.index] = grp + "_s" + d["subcluster"].astype(str)
    A.obs["sub"] = subcol.reindex(A.obs_names).values
    sym2gid = {}
    for gid, sym in zip(A.var_names.astype(str), A.var["sym"].astype(str)):
        sym2gid.setdefault(sym, gid)
    for name, genes in PROGRAMS.items():
        ids = [sym2gid[g] for g in genes if g in sym2gid]
        sc.tl.score_genes(A, ids, score_name=name, random_state=0)
    rib = [gid for gid, s in zip(A.var_names.astype(str), A.var["sym"].astype(str))
           if s.lower().startswith(("rpl", "rps"))]
    sc.tl.score_genes(A, rib, score_name="Ribosome", random_state=0)
    score_cols = list(PROGRAMS.keys()) + ["Ribosome"]

    # top markers
    top = {}
    for f in glob.glob(os.path.join(SUB_DIR, "*_subcluster_markers.csv")):
        d = pd.read_csv(f)
        for cl, sub in d.groupby("subcluster"):
            top[f"{os.path.basename(f).split('_')[0]}_s{cl}"] = ", ".join(sub.head(5).gene.astype(str))

    rows = []
    for r in comp.itertuples():
        key = r.key; parent = str(r.population)
        nctrl, nhsp = int(r.n_control), int(r.n_hemospasia)
        total = nctrl + nhsp
        parent_total = int(comp[comp.population == parent][["n_control", "n_hemospasia"]].sum().sum())
        parent_min = int(comp[comp.population == parent][["n_control", "n_hemospasia"]].sum().min())
        if parent in ("MyoEpi", "M") or parent_total < 100 or parent_min < 30:
            status = "excluded_parent_small"
        elif total < 20 or min(nctrl, nhsp) < 5:
            status = "excluded_subcluster_small"
        else:
            status = "eligible"
        sub = A[A.obs["sub"].astype(str) == key]
        score = {c: float(sub.obs[c].mean()) for c in score_cols} if len(sub) else {c: np.nan for c in score_cols}
        fold = float(r.fold_hemospasia_vs_control) if pd.notna(r.fold_hemospasia_vs_control) else np.nan
        if status == "eligible":
            change = "expanded" if fold >= 1.25 else ("reduced" if fold <= 0.75 else "stable")
        else:
            change = "not_analyzed"
        rows.append({
            "parent": parent, "subcluster": str(r.subcluster), "key": key,
            "proposed_state_label": LABELS.get(key, "unlabeled"),
            "n_control": nctrl, "n_hemospasia": nhsp, "n_total": total,
            "pct_control": float(r.pct_control), "pct_hemospasia": float(r.pct_hemospasia),
            "fold_hemospasia_vs_control": fold, "change": change,
            "analysis_status": status, "top_markers": top.get(key, ""),
            **score,
        })
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "Fig5_dynamic_subcluster_annotation_table.csv"),
             index=False, encoding="utf-8-sig")

    E = T[T.analysis_status == "eligible"].copy()
    E["order"] = E.key.map({k: i for i, k in enumerate(LABEL_ORDER)})
    E = E.sort_values("order").reset_index(drop=True)

    # Annotation heatmap
    heat_cols = ["PH/progenitor", "Macrophage/immune", "Absorptive/metabolic",
                 "ECM/remodeling", "Stress", "PM-like/cuticle", "Ribosome", "G2M"]
    H = E[heat_cols].astype(float)
    Z = (H - H.mean(axis=0)) / H.std(axis=0, ddof=0).replace(0, 1)
    fig = plt.figure(figsize=(7.2, 4.6), facecolor="white")
    ax = fig.add_axes([0.43, 0.12, 0.40, 0.76])
    im = ax.imshow(np.clip(Z.values, -2, 2), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(heat_cols)))
    ax.set_xticklabels(heat_cols, rotation=42, ha="right", rotation_mode="anchor", fontsize=5.2)
    ax.set_yticks(np.arange(len(E)))
    ax.set_yticklabels([f"{k}  |  {lab}" for k, lab in zip(E.key, E.proposed_state_label)], fontsize=5.0)
    ax.tick_params(axis="both", length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # axes-level separators
    for i in range(1, len(E)):
        if E.loc[i, "parent"] != E.loc[i - 1, "parent"]:
            ax.axhline(i - 0.5, color="#444444", lw=0.7)
    cax = fig.add_axes([0.85, 0.30, 0.018, 0.38])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Program score (row z-score)", fontsize=5.2, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.3, width=0.4)
    cb.outline.set_linewidth(0.4)
    fig.text(0.50, 0.96, "Functional programs of eligible dynamic substates", ha="center",
             va="top", fontsize=7.2)
    fig.text(0.50, 0.03, "Only substates with >=20 cells and >=5 cells per group are shown. "
             "Labels are exploratory and based on marker and program-score evidence.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out1 = os.path.join(PARTS, "Fig5_dynamic_subcluster_annotation_heatmap_preview.png")
    fig.savefig(out1, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out1)
    plt.close(fig)

    # Paired composition dot plot
    fig = plt.figure(figsize=(7.2, 4.6), facecolor="white")
    ax = fig.add_axes([0.43, 0.12, 0.46, 0.76])
    y = np.arange(len(E))[::-1]
    for yy, r in zip(y, E.itertuples()):
        x1, x2 = float(r.pct_control), float(r.pct_hemospasia)
        ax.plot([x1, x2], [yy, yy], color="#BDBDBD", lw=0.9, zorder=1)
        ax.scatter([x1], [yy], s=18, c="#4C78A8", edgecolors="white", linewidths=0.35, zorder=3)
        ax.scatter([x2], [yy], s=18, c="#D6604D", edgecolors="white", linewidths=0.35, zorder=3)
        ax.text(104, yy, f"{r.fold_hemospasia_vs_control:.2f}×", ha="left", va="center",
                fontsize=4.9, color="#B2182B" if r.fold_hemospasia_vs_control > 1 else "#2166AC")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k}  |  {lab}" for k, lab in zip(E.key, E.proposed_state_label)], fontsize=5.0)
    ax.set_xlim(0, 118); ax.set_ylim(-0.6, len(E) - 0.4)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Percentage within group (%)", fontsize=5.5)
    ax.tick_params(axis="x", labelsize=5.0, length=1.4, width=0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#EEEEEE", lw=0.35, zorder=0)
    for i in range(1, len(E)):
        if E.loc[i, "parent"] != E.loc[i - 1, "parent"]:
            ax.axhline(y[i] + 0.5, color="#444444", lw=0.7)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="#4C78A8", markersize=4, label="Control"),
                       Line2D([0], [0], marker="o", color="none", markerfacecolor="#D6604D", markersize=4, label="Hemospasia")],
              loc="upper right", bbox_to_anchor=(1.0, 1.12), ncol=2, frameon=False, fontsize=5.0)
    fig.text(0.06, 0.965, "Subcluster composition change", ha="left", va="top", fontsize=7.2)
    fig.text(0.50, 0.03, "Fold change is shown at the right. Only eligible substates are shown.",
             ha="center", va="bottom", fontsize=5.0, color="#555555")
    out2 = os.path.join(PARTS, "Fig5_dynamic_subcluster_change_dotplot_preview.png")
    fig.savefig(out2, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out2)
    plt.close(fig)


if __name__ == "__main__":
    main()
