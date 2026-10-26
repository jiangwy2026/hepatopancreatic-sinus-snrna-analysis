# -*- coding: utf-8 -*-
"""
146_score_published_hep_markers.py
用 PMID 42004253 (PMC13089082) 公布的 12 类肝胰腺细胞 marker，对我们 10 个群打分 -> 标准命名对照
"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
CL = ["0","1","2","3","4","5","6","7","8","9"]
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NEU","9":"NEU"}
SETS = {
 "F cell": ["Amy2","Apod","CTSC","LCP2"],
 "R cell": ["Slc37a2","ACAC","SCD5"],
 "B cell": ["Gs2","vha-10","LOC113800236"],
 "E cell": ["Dsp1","eIF-5A","EIF5A"],
 "M cell": ["ChAT","ADCY5"],
 "endothelial": ["Igfbp7","SPON1"],
 "granulocyte": ["PEN-3","LOC113801825"],
 "phagocyte": ["PRDX5","Nlrp3","ALF"],
 "immature R": ["NAAT1","Picot"],
 "keratinocyte": ["Clca4a","Cd63"],
 "myoepithelial": ["Plg","LanB2"],
 "peritrophic matrix": ["MCFD2","FSTL4"],
}
def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sym = (A.var["sym"] if "sym" in A.var.columns else A.var["symbol"]).astype(str).values; ids = A.var_names.astype(str).values
    tot = np.asarray(A.layers["counts"].sum(0)).ravel()
    best = {}
    for i, s in enumerate(sym):
        if s not in best or tot[i] > tot[best[s]]: best[s] = i
    cl = A.obs["vendor_cluster"].astype(str).values
    own = A.obs_names

    print("=== marker 基因可用性 ===")
    used = {}; missing_total = []
    for name, gl in SETS.items():
        u = [g for g in gl if g in best]
        if name == "E cell" and "Dsp1" in u and ("eIF-5A" in u or "EIF5A" in u):
            u = [g for g in u if g != ("EIF5A" if "eIF-5A" in u else "eIF-5A")]
        used[name] = u; miss = [g for g in gl if g not in best]
        missing_total += miss
        print("  %-20s 使用 %-2d 个: %-40s 缺失: %s" % (name, len(u), ", ".join(u), ", ".join(miss) if miss else "-"))
    print("  缺失基因共 %d 个: %s" % (len(set(missing_total)), sorted(set(missing_total))))

    M = pd.DataFrame(index=list(SETS), columns=CL, dtype=float)
    for name, u in used.items():
        gi = [best[g] for g in u]
        if not gi: continue
        for c in CL:
            m = cl == c
            M.loc[name, c] = float(np.asarray(A[m][:, gi].X.mean()))
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
    print("\n=== 各群 × 12 类 published marker 的 z-score ===")
    show = Z.copy(); show.columns = ["%s(C%s)" % (OUR[c], c) for c in CL]
    print(show.round(2).to_string())
    Z.to_csv(os.path.join(OUT, "hep12_marker_scores_z.csv"))
    M.to_csv(os.path.join(OUT, "hep12_marker_scores_raw.csv"))

    print("\n=== 每个群的 top-3 匹配 ===")
    rows = []
    for c in CL:
        col = Z[c].dropna().sort_values(ascending=False)
        rows.append({"cluster": "C"+c, "our_name": OUR[c],
                     "top1": "%s (%.2f)" % (col.index[0], col.iloc[0]),
                     "top2": "%s (%.2f)" % (col.index[1], col.iloc[1]),
                     "top3": "%s (%.2f)" % (col.index[2], col.iloc[2])})
        print("  C%-2s %-4s -> %s | %s | %s" % (c, OUR[c], rows[-1]["top1"], rows[-1]["top2"], rows[-1]["top3"]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "hep12_top_match_per_cluster.csv"), index=False)

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    vmax = np.nanmax(np.abs(Z.values))
    im = ax.imshow(Z.values, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax), aspect="auto")
    ax.set_xticks(range(len(CL))); ax.set_xticklabels(["%s\n(C%s)" % (OUR[c], c) for c in CL], fontsize=6)
    ax.set_yticks(range(len(Z))); ax.set_yticklabels(Z.index, fontsize=6.5)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            v = Z.values[i, j]
            if np.isfinite(v) and abs(v) > 1.2:
                ax.text(j, i, "%.1f" % v, fontsize=5, ha="center", va="center", color="white" if abs(v) > 2 else "#333")
    ax.set_title("Published hepatopancreas cell-type markers\n(PMID 42004253, 12 cell types) vs our clusters", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("z-score", fontsize=6)
    cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
    for s in ax.spines.values(): s.set_visible(False)
    p = os.path.join(FIG, "Fig_hep12_marker_mapping_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("\nsaved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\hep12_mapping.log", "w", encoding="utf8", buffering=1)
    main()
