# -*- coding: utf-8 -*-
"""
150_bubble_annotation.py
按用户设想：气泡大小 = 在该群中被"支持"的 marker 基因数；颜色 = 匹配准确度（one-vs-rest AUROC 中位数）
左：eLife 血细胞 6 亚型；右：肝胰腺 12 类
"""
import io, sys, os, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.stats import rankdata
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
CL = ["0","1","2","3","4","5","6","7","8","9"]
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"Neu8","9":"Neu9"}
EMAP = {"elife-80127-supp3-v2.csv":"PH1","elife-80127-supp4-v2.csv":"GH2",
        "elife-80127-supp5-v2.csv":"MH2","elife-80127-supp6-v2.csv":"PH2",
        "elife-80127-supp7-v2.csv":"GH1","elife-80127-supp8-v2.csv":"MH1"}
HP = {"F": ["Amy2","Apod","CTSC","LCP2"], "R": ["Slc37a2","ACAC","SCD5"],
      "B": ["Gs2","vha-10","LOC113800236"], "E": ["Dsp1","eIF-5A"],
      "M": ["ChAT","ADCY5"], "endothelial": ["Igfbp7","SPON1"],
      "granulocyte": ["PEN-3","LOC113801825"], "phagocyte": ["PRDX5","Nlrp3","ALF"],
      "immatureR": ["NAAT1","Picot"], "keratinocyte": ["Clca4a","Cd63"],
      "myoepithelial": ["Plg","LanB2"], "peritrophic": ["MCFD2","FSTL4"]}

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    ids = list(A.var_names.astype(str)); id2i = {g: i for i, g in enumerate(ids)}
    sym = A.var["sym"].astype(str).values
    tg = np.asarray(A.layers["counts"].sum(0)).ravel()
    best = {}
    for i, s in enumerate(sym):
        if s not in best or tg[i] > tg[best[s]]: best[s] = i
    cl = A.obs["vendor_cluster"].astype(str).values
    n = A.n_obs

    # 参考集合
    sets = {}
    for f, name in EMAP.items():
        d = pd.read_csv(os.path.join(OUT, "elife_ref", f))
        if "Log2FC" in d.columns: d = d.sort_values("Log2FC", ascending=False)
        sets["eLife " + name] = [g for g in d["Gene ID"].astype(str) if g in id2i][:100]
    for k, gl in HP.items():
        sets["HP " + k] = [ids[best[g]] for g in gl if g in best]
    allg = sorted({g for v in sets.values() for g in v})
    gi = [id2i[g] for g in allg]
    X = A[:, gi].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    # pct 矩阵 + 每个基因的 rank
    pct = np.zeros((len(allg), len(CL)))
    R = np.zeros_like(X)
    for j, c in enumerate(CL):
        m = cl == c
        pct[:, j] = 100 * (X[m] > 0).mean(0)
    for k in range(X.shape[1]):
        R[:, k] = rankdata(X[:, k])
    gpos = {g: k for k, g in enumerate(allg)}
    # AUROC per gene per cluster
    AUC = np.zeros((len(allg), len(CL)))
    for j, c in enumerate(CL):
        m = cl == c; npos = m.sum(); nneg = n - npos
        s = R[m].sum(0)
        AUC[:, j] = (s - npos * (npos + 1) / 2) / (npos * nneg)

    rows = []
    out = {}
    for name, gl in sets.items():
        kk = [gpos[g] for g in gl if g in gpos]
        for j, c in enumerate(CL):
            p_c = pct[kk, j]; p_o = np.max(np.delete(pct[kk], j, axis=1), axis=1)
            sup = (p_c > p_o) & (p_c > 5)
            out[(name, c)] = (int(sup.sum()), float(np.median(AUC[kk, j])))
    names = list(sets)
    tab_n = pd.DataFrame([[out[(nm, c)][0] for nm in names] for c in CL], index=["%s(C%s)" % (OUR[c], c) for c in CL], columns=names)
    tab_a = pd.DataFrame([[out[(nm, c)][1] for nm in names] for c in CL], index=tab_n.index, columns=names)
    tab_n.to_csv(os.path.join(OUT, "bubble_n_supported.csv")); tab_a.to_csv(os.path.join(OUT, "bubble_auc.csv"))
    print("=== 支持基因数 ==="); print(tab_n.to_string())
    print("\n=== 准确度 (median AUROC) ==="); print(tab_a.round(3).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.9))
    for ax, cols, title in ((axes[0], [n_ for n_ in names if n_.startswith("eLife")], "eLife hemocyte subtypes"),
                            (axes[1], [n_ for n_ in names if n_.startswith("HP")], "Hepatopancreas cell types")):
        N = tab_n[cols].values; A_ = tab_a[cols].values
        for i in range(N.shape[0]):
            for j in range(N.shape[1]):
                if N[i, j] == 0: continue
                ax.scatter(j, i, s=8 + N[i, j] * 3.2, c=[plt.cm.viridis((np.clip(A_[i, j], 0.5, 1.0) - 0.5) / 0.5)],
                           edgecolors="#555", linewidths=0.3, zorder=3)
        ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(len(CL) - 0.4, -0.6)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace("eLife ", "").replace("HP ", "") for c in cols],
                                                            rotation=45, ha="right", fontsize=6)
        ax.set_yticks(range(len(CL))); ax.set_yticklabels(tab_n.index, fontsize=6.5)
        ax.set_title(title, fontsize=8)
        ax.grid(True, ls=(0, (1, 2)), lw=0.3, color="#E5E5E5", zorder=0)
        for s in ax.spines.values(): s.set_visible(False)
    axes[1].set_yticklabels([])
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=Normalize(0.5, 1.0))
    cb = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
    cb.set_label("Matching accuracy (median AUROC)", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)
    for s, nn in ((8, 2), (16, 8), (24, 20)):
        fig.text(0.012, 0.62 - s * 0.006, "●", fontsize=3 + s * 0.28, color="#888")
    fig.suptitle("Annotation evidence: bubble size = supported marker genes, colour = median one-vs-rest AUROC",
                 fontsize=8.5, y=1.0)
    p = os.path.join(FIG, "Fig_annotation_bubbles_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("\nsaved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\annotation_bubbles.log", "w", encoding="utf8", buffering=1)
    main()
