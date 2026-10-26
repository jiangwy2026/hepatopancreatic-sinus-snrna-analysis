# -*- coding: utf-8 -*-
"""
148_specificity_figure.py
更准确的注释证据：(A) 各参考细胞类型的"marker 特异性"（中位检出率差）与 marker 数
(B) 关键基因逐群检出率点图 —— 直接看出高分是被 1-2 个广谱基因撑起来的
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
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"neu8","9":"neu9"}
HP = {"F cell": ["Amy2","Apod","CTSC","LCP2"], "R cell": ["Slc37a2","ACAC","SCD5"],
      "B cell": ["Gs2","vha-10","LOC113800236"], "E cell": ["Dsp1","eIF-5A"],
      "M cell": ["ChAT","ADCY5"], "endothelial": ["Igfbp7","SPON1"],
      "granulocyte": ["PEN-3","LOC113801825"], "phagocyte": ["PRDX5","Nlrp3","ALF"],
      "immature R": ["NAAT1","Picot"], "keratinocyte": ["Clca4a","Cd63"],
      "myoepithelial": ["Plg","LanB2"], "peritrophic matrix": ["MCFD2","FSTL4"]}

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sym = A.var["sym"].astype(str).values; ids = A.var_names.astype(str).values
    tg = np.asarray(A.layers["counts"].sum(0)).ravel()
    best = {}
    for i, s in enumerate(sym):
        if s not in best or tg[i] > tot if False else False: pass
    best = {}
    for i, s in enumerate(sym):
        if s not in best or tg[i] > tg[best[s]]: best[s] = i
    id2i = {g: i for i, g in enumerate(ids)}
    cl = A.obs["vendor_cluster"].astype(str).values

    def profiles(gene_idx):
        pct = np.zeros((len(gene_idx), len(CL))); mean = np.zeros_like(pct)
        for j, c in enumerate(CL):
            m = cl == c
            X = A[m][:, gene_idx].X
            X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
            pct[:, j] = 100 * (X > 0).mean(0); mean[:, j] = X.mean(0)
        return pct, mean

    # eLife marker sets
    elife = {}
    for f in sorted(glob.glob(os.path.join(OUT, "elife_ref", "elife-80127-supp*-v2.csv"))):
        d = pd.read_csv(f)
        if "Target_Cluster" not in d.columns: continue
        tc = str(d["Target_Cluster"].iloc[0])
        if tc not in ("PH1","PH2","GH1","GH2","MH1","MH2"): continue
        d = d[(d.get("Qvalue", 0) < 0.05)] if "Qvalue" in d.columns else d
        if "Log2FC" in d.columns: d = d.sort_values("Log2FC", ascending=False)
        gids = [g for g in d["Gene ID"].astype(str) if g in id2i][:120]
        elife[tc] = gids
    print("eLife marker 数:", {k: len(v) for k, v in elife.items()})

    rows = []
    def add_set(label, lineage, gids):
        if not gids: return None
        pct, mean = profiles([id2i[g] for g in gids])
        mc = int(np.argmax(mean.mean(axis=0)))
        gap = pct[:, mc] - np.max(np.delete(pct, mc, axis=1), axis=1)
        rows.append({"reference": label, "lineage": lineage, "n_markers": len(gids),
                     "matched_cluster": "%s(C%s)" % (OUR[CL[mc]], CL[mc]),
                     "median_gap_pp": float(np.median(gap)), "frac_gap_gt20": float((gap > 20).mean()),
                     "median_pct_matched": float(np.median(pct[:, mc]))})
        return pct, mean, gap, mc
    for k, v in elife.items():
        add_set("eLife " + k, "blood", v)
    for k, gl in HP.items():
        gids = [ids[best[g]] for g in gl if g in best]
        add_set("HP " + k, "hepato", gids)
    R = pd.DataFrame(rows).sort_values("median_gap_pp", ascending=False)
    print("\n=== 各参考细胞类型的 marker 特异性 ===")
    print(R.round(1).to_string(index=False))
    R.to_csv(os.path.join(OUT, "annotation_specificity_table.csv"), index=False)

    # ---- 图：左 = 特异性条形图；右 = 关键基因点图
    genes = [("Igfbp7","endothelial"),("SPON1","endothelial"),("Dsp1","E cell"),("eIF-5A","E cell"),
             ("Clca4a","keratinocyte"),("Cd63","keratinocyte"),("LOC113801825","granulocyte"),
             ("PEN-3","granulocyte"),("Nlrp3","phagocyte"),("ALF","phagocyte"),
             ("ncbi_113814652","PH (ccnd2)"),("ncbi_113823802","PH (Grik2)"),
             ("ncbi_113823945","PH (TGM)"),("ncbi_113812814","PH (PPAF2)")]
    gidx, glab = [], []
    for g, lab in genes:
        if g in id2i: gidx.append(id2i[g]); glab.append("%s [%s]" % (sym[id2i[g]], lab))
        elif g in best: gidx.append(best[g]); glab.append("%s [%s]" % (g, lab))
    pct, mean = profiles(gidx)
    z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)

    fig = plt.figure(figsize=(11.2, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.42, left=0.03, right=0.97, top=0.90, bottom=0.20)
    ax = fig.add_subplot(gs[0])
    colors = ["#B2182B" if l == "blood" else "#2166AC" for l in R.lineage]
    y = np.arange(len(R))
    ax.barh(y, R.median_gap_pp, color=colors, height=0.68)
    ax.set_yticks(y); ax.set_yticklabels(["%s (n=%d)" % (r.reference, r.n_markers) for r in R.itertuples()], fontsize=6)
    ax.invert_yaxis(); ax.axvline(0, color="#888", lw=0.6)
    ax.set_xlabel("Marker specificity: median(% in matched cluster − max % elsewhere), pp", fontsize=6.5)
    ax.set_title("How specific are the markers behind each match?", fontsize=8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=5.5)
    ax = fig.add_subplot(gs[1])
    for i in range(len(glab)):
        for j in range(len(CL)):
            if pct[i, j] <= 1: continue
            ax.scatter(j, i, s=6 + pct[i, j] / 100 * 80, c=[plt.cm.RdBu_r((np.clip(z[i, j], -2, 2) + 2) / 4)],
                       edgecolors="#555", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.6, len(CL) - 0.4); ax.set_ylim(len(glab) - 0.4, -0.6)
    ax.set_xticks(range(len(CL))); ax.set_xticklabels(["%s(C%s)" % (OUR[c], c) for c in CL], rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(glab))); ax.set_yticklabels(glab, fontsize=5.6)
    ax.grid(True, ls=(0, (1, 2)), lw=0.3, color="#DDD", zorder=0)
    ax.set_title("Gene-level evidence (dot = % expressing, colour = scaled mean)", fontsize=8)
    for s in ax.spines.values(): s.set_visible(False)
    fig.text(0.03, 0.95, "red = hemocyte reference   blue = hepatopancreas reference", fontsize=6.5, color="#555")
    p = os.path.join(FIG, "Fig_annotation_specificity_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("\nsaved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\annotation_specificity.log", "w", encoding="utf8", buffering=1)
    main()
