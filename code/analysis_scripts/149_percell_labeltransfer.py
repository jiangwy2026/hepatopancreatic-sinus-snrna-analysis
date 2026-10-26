# -*- coding: utf-8 -*-
"""
149_percell_labeltransfer.py
逐细胞竞争性标签转移（score_genes + 表达匹配对照基因）-> 混淆矩阵
比"群均值 z-score"更稳健：每个细胞只能选一个身份，广谱基因无法靠数量取胜
"""
import io, sys, os, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
CL = ["0","1","2","3","4","5","6","7","8","9"]
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"C8","9":"C9"}
BLOOD_CL = ["0","1","2","3"]; HP_CL = ["4","5","6","7","8","9"]
EMAP = {"elife-80127-supp3-v2.csv":"PH1","elife-80127-supp4-v2.csv":"GH2",
        "elife-80127-supp5-v2.csv":"MH2","elife-80127-supp6-v2.csv":"PH2",
        "elife-80127-supp7-v2.csv":"GH1","elife-80127-supp8-v2.csv":"MH1"}
HP = {"F cell": ["Amy2","Apod","CTSC","LCP2"], "R cell": ["Slc37a2","ACAC","SCD5"],
      "B cell": ["Gs2","vha-10","LOC113800236"], "E cell": ["Dsp1","eIF-5A"],
      "M cell": ["ChAT","ADCY5"], "endothelial": ["Igfbp7","SPON1"],
      "granulocyte": ["PEN-3","LOC113801825"], "phagocyte": ["PRDX5","Nlrp3","ALF"],
      "immature R": ["NAAT1","Picot"], "keratinocyte": ["Clca4a","Cd63"],
      "myoepithelial": ["Plg","LanB2"], "peritrophic matrix": ["MCFD2","FSTL4"]}

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

    # ---- eLife 6 亚型打分
    sets = {}
    for f, name in EMAP.items():
        d = pd.read_csv(os.path.join(OUT, "elife_ref", f))
        if "Log2FC" in d.columns: d = d.sort_values("Log2FC", ascending=False)
        sets[name] = [g for g in d["Gene ID"].astype(str) if g in id2i][:100]
    print("eLife 打分基因数:", {k: len(v) for k, v in sets.items()})
    S = pd.DataFrame(index=A.obs_names)
    for k, gl in sets.items():
        sc.tl.score_genes(A, gl, ctrl_size=100, score_name="s_" + k, random_state=0)
        S[k] = A.obs["s_" + k].values
    S.to_csv(os.path.join(OUT, "elife_percell_scores.csv"))
    asg = S.idxmax(1); conf = (S.max(1) - S.apply(lambda r: r.nlargest(2).iloc[-1], axis=1))
    tab = pd.crosstab(pd.Series(cl, index=A.obs_names, name="our"), asg, normalize="index") * 100
    tab = tab.reindex([c for c in CL if c in tab.index])
    print("\n=== eLife 6 亚型：逐细胞标签转移（行 = 我们的群，列 = 被分配的 eLife 亚型，%）===")
    print(tab.round(1).to_string())

    # ---- 肝胰腺 12 类打分
    HPs = {k: [ids[best[g]] for g in gl if g in best] for k, gl in HP.items()}
    S2 = pd.DataFrame(index=A.obs_names)
    for k, gl in HPs.items():
        if not gl: continue
        sc.tl.score_genes(A, gl, ctrl_size=100, score_name="h_" + k.replace(" ", "_"), random_state=0)
        S2[k] = A.obs["h_" + k.replace(" ", "_")].values
    asg2 = S2.idxmax(1)
    tab2 = pd.crosstab(pd.Series(cl, index=A.obs_names, name="our"), asg2, normalize="index") * 100
    tab2 = tab2.reindex([c for c in HP_CL if c in tab2.index])
    print("\n=== 肝胰腺 12 类：逐细胞标签转移（%）===")
    print(tab2.round(1).to_string())

    # ---- 图
    fig = plt.figure(figsize=(11.0, 3.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.42, left=0.07, right=0.97, top=0.86, bottom=0.26)
    ax = fig.add_subplot(gs[0])
    M1 = tab.values; im1 = ax.imshow(M1, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(tab.shape[1])); ax.set_xticklabels(tab.columns, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(tab.shape[0])); ax.set_yticklabels(["%s(C%s)" % (OUR[c], c) for c in tab.index], fontsize=7)
    for i in range(M1.shape[0]):
        for j in range(M1.shape[1]):
            if M1[i, j] >= 5: ax.text(j, i, "%.0f" % M1[i, j], fontsize=5.5, ha="center", va="center",
                                      color="white" if M1[i, j] > 55 else "#333")
    ax.set_title("Per-cell assignment to eLife hemocyte subtypes\n(row-normalised %)", fontsize=7.5)
    ax.set_xlabel("Assigned eLife subtype", fontsize=6.5)
    ax = fig.add_subplot(gs[1])
    M2 = tab2.values; im2 = ax.imshow(M2, cmap="Greens", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(tab2.shape[1])); ax.set_xticklabels(tab2.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(tab2.shape[0])); ax.set_yticklabels(["%s(C%s)" % (OUR[c], c) for c in tab2.index], fontsize=7)
    for i in range(M2.shape[0]):
        for j in range(M2.shape[1]):
            if M2[i, j] >= 5: ax.text(j, i, "%.0f" % M2[i, j], fontsize=5.2, ha="center", va="center",
                                      color="white" if M2[i, j] > 55 else "#333")
    ax.set_title("Per-cell assignment to hepatopancreas cell types\n(row-normalised %)", fontsize=7.5)
    ax.set_xlabel("Assigned cell type", fontsize=6.5)
    p = os.path.join(FIG, "Fig_annotation_confusion_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("\nsaved", p)
    tab.to_csv(os.path.join(OUT, "confusion_elife_pct.csv")); tab2.to_csv(os.path.join(OUT, "confusion_hepato_pct.csv"))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\confusion.log", "w", encoding="utf8", buffering=1)
    main()
