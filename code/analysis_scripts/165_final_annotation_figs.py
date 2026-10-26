# -*- coding: utf-8 -*-
"""
165_final_annotation_figs.py
1) 统一新命名：C0=TS, C4=PM-like, C8=MyoEpi, C9=M；C1=GH, C2=PH, C3=MH2, C5=F, C6=R, C7=B
2) 用新名字重画两张气泡图（色标缩小、肝胰腺面板气泡放大）
3) 保存程序表达均值表（新命名）
4) 输出备选参考文献清单
"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")

OLD = ["MH1(C0)","GH(C1)","PH(C2)","MH2(C3)","Epi(C4)","F(C5)","R(C6)","B(C7)","Neu8(C8)","Neu9(C9)"]
NEW = {"MH1(C0)":"TS(C0)","GH(C1)":"GH(C1)","PH(C2)":"PH(C2)","MH2(C3)":"MH2(C3)",
       "Epi(C4)":"PM-like(C4)","F(C5)":"F(C5)","R(C6)":"R(C6)","B(C7)":"B(C7)",
       "Neu8(C8)":"MyoEpi(C8)","Neu9(C9)":"M(C9)"}
ORDER = ["PH(C2)","GH(C1)","TS(C0)","MH2(C3)","F(C5)","R(C6)","B(C7)","M(C9)","MyoEpi(C8)","PM-like(C4)"]

def draw(n_tab, a_tab, cols, title, fname, figsize, scale, legend, xfs=6.0, rot=45, lx=0.45):
    n = n_tab.loc[ORDER, cols].values; a = a_tab.loc[ORDER, cols].values
    fig, ax = plt.subplots(figsize=figsize)
    for i in range(n.shape[0]):
        for j in range(n.shape[1]):
            if not np.isfinite(a[i, j]): continue
            ax.scatter(j, i, s=14 + n[i, j] * scale,
                       c=[plt.cm.viridis((np.clip(a[i, j], 0.5, 1.0) - 0.5) / 0.5)],
                       edgecolors="#555555", linewidths=0.3, zorder=3)
    ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(len(ORDER) - 0.4, -0.6)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("eLife ", "").replace("HP ", "") for c in cols], rotation=rot, ha="right", fontsize=xfs)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER, fontsize=7)
    ax.set_title(title, fontsize=8.5, pad=6)
    ax.grid(True, ls=(0, (1, 2)), lw=0.3, color="#E5E5E5", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=Normalize(0.5, 1.0))
    cb = fig.colorbar(sm, ax=ax, fraction=0.022, shrink=0.58, pad=0.02)
    cb.set_label("AUROC", fontsize=5.5); cb.ax.tick_params(labelsize=4.5, length=1.2)
    cb.outline.set_visible(False)
    ax.text(len(cols) + lx + 0.02, -0.30, "genes", fontsize=5.8, va="center", clip_on=False)
    for k, v in enumerate(legend):
        ax.scatter(len(cols) + lx, 0.35 + k * 1.05, s=14 + v * scale, c="none", edgecolors="#666666",
                   linewidths=0.4, clip_on=False)
        ax.text(len(cols) + lx + 0.32, 0.35 + k * 1.05, "%d" % v, fontsize=5.8, va="center", clip_on=False)
    p = os.path.join(FIG, fname)
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

def main():
    N = pd.read_csv(os.path.join(OUT, "bubble2_n_supported.csv"), index_col=0).rename(index=NEW)
    A = pd.read_csv(os.path.join(OUT, "bubble2_auc.csv"), index_col=0).rename(index=NEW)
    els = [c for c in N.columns if c.startswith("eLife")]
    hps = [c for c in N.columns if c.startswith("HP")]
    draw(N, A, els, "eLife hemocyte subtypes", "Fig_annotation_bubbles_blood_preview.png",
         (4.8, 3.4), 3.6, (5, 20, 50))
    draw(N, A, hps, "Hepatopancreas cell types (PMID 42004253)", "Fig_annotation_bubbles_hepatopancreas_preview.png",
         (6.9, 3.4), 30, (1, 2, 4), xfs=5.8, lx=1.15)

    # 程序表（新命名）
    P = pd.read_csv(os.path.join(OUT, "C0_metabolic_programs.csv"), index_col=0)
    PMAP = {"MH1(C0)": "TS", "GH": "GH", "PH": "PH", "MH2": "MH2", "Epi": "PM-like",
            "F": "F", "R": "R", "B": "B", "Neu8": "MyoEpi", "Neu9": "M"}
    P = P.rename(columns=PMAP)
    P = P[["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]]
    P.index.name = "program"
    P.to_csv(os.path.join(OUT, "metabolic_programs_by_cluster.csv"))
    # 热图（行内 z-score）
    Z = P.sub(P.mean(1), axis=0).div(P.std(1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(5.6, 2.2))
    im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(P.shape[1])); ax.set_xticklabels(P.columns, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(P.shape[0])); ax.set_yticklabels(P.index, fontsize=6.2)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if abs(Z.values[i, j]) > 1.0:
                ax.text(j, i, "%.1f" % Z.values[i, j], fontsize=5, ha="center", va="center",
                        color="white" if abs(Z.values[i, j]) > 1.6 else "#333")
    ax.set_title("Metabolic / biosynthetic programs by cluster (row z-score)", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, shrink=0.7, pad=0.02)
    cb.ax.tick_params(labelsize=4.8); cb.outline.set_visible(False)
    for s in ax.spines.values(): s.set_visible(False)
    fig.savefig(os.path.join(FIG, "Fig_metabolic_programs_preview.png"), dpi=400, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("\n=== 程序表达均值表（新命名）===")
    print(P.to_string())

    refs = """# 亚群命名相关参考文献（候选）

## 已用于命名（主要依据）

1. **Yang et al., 2022, eLife 11:e80127** — *Single-cell RNA sequencing analysis of shrimp immune cells identifies macrophage-like phagocytes*（血细胞 6 亚型 PH1/PH2/GH1/GH2/MH1/MH2 的 marker 表）
2. **PMID 42004253**, *Anim Nutr* 2026;25:458–475（PMC13089082）— 肝胰腺 12 类细胞 marker（F/R/B/E/M/endothelial/granulocyte/phagocyte/immature R/keratinocyte/myoepithelial/peritrophic matrix）
3. **PMID 42000089**, *Fish Shellfish Immunol* 2026;174:111364 — 肝胰腺 5 类细胞（B/F/E/M/R）
4. **GSE306674 / PMID 41694335**, *Front Immunol* 2026;17:1713369 — 肝胰腺单细胞图谱（本项目用于标签转移的参考数据）

## 备选/旁证（数据未公开或未使用其 marker）

5. **PMID 42250742**, *Fish Shellfish Immunol* 2026;176:111494 — 血细胞图谱，描述"high-metabolic biosynthetic transition state"（用于支持 C0 → **TS** 的命名；其数据与补充材料未公开）
6. **PMID 41638371** — 香蕉虾血细胞 8 群，YHV 感染后未成熟祖细胞群扩增（emergency hematopoiesis）
7. **PMID 41248710**, *J Adv Res* 2026 — 淋巴器官 Oka 单核图谱（8 类：progenitors / macrophage-like / lymphocyte precursor-like 等）
8. **PMID 42044139**, *PLoS One* 2026 — 首个对虾淋巴器官单核细胞图谱
9. **PMID 38323810**, *J Virol* 2024 — 凡纳滨对虾血细胞 WSSV 单细胞
10. **PMID 36761772**, *Front Immunol* 2023 — 甲壳类血细胞 11/12 类
11. **PMID 41038398**, *Genomics* 2025 — 日本对虾血细胞 13 亚群

## 命名方案（本数据最终）

| 群 | 名称 | 依据 |
|---|---|---|
| C2 | **PH** (prohemocyte) | eLife PH1/PH2（AUROC 0.89/0.87） |
| C1 | **GH** (granular hemocyte) | eLife GH1/GH2（0.85/0.84） |
| C0 | **TS** (biosynthetic transition state) | 核糖体/翻译/OXPHOS 在血细胞中最高；与 PMID 42250742 描述一致 |
| C3 | **MH2** (macrophage-like phagocyte) | eLife MH2（0.96）+ phagocyte markers |
| C5 | **F** | HP F cell（0.84） |
| C6 | **R** | HP R cell（0.64） |
| C7 | **B** | HP B cell（0.73） |
| C9 | **M** | HP M cell（0.73） |
| C8 | **MyoEpi** (myoepithelial-like) | HP myoepithelial（0.59，弱） |
| C4 | **PM-like** (peritrophic matrix-like) | HP peritrophic（0.53）与 R（0.62）之间，身份未定 |
"""
    open(os.path.join(OUT, "REFERENCES_ANNOTATION.md"), "w", encoding="utf8").write(refs)
    print("\nsaved REFERENCES_ANNOTATION.md")

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\final_annotation.log", "w", encoding="utf8", buffering=1)
    main()
