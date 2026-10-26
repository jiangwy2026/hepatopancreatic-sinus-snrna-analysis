# -*- coding: utf-8 -*-
"""250_step3_lyo_program_mirror.py

Step 3: lymphoid-organ top50 programs scored in the hepatic-sinusoid/hemocyte
dataset with a UCell-style rank score, plus hypergeometric overlap tests.

Inputs:
  DeGene.list.xls (Top50.DeGene.list, 16 LyO clusters)
  step2/step2_integrated.h5ad (9,658 cells, lognorm)
  markers_10groups.csv (top markers per cell type)
Outputs in analysis_output/step3/.
"""
import os, json, sys
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\桌面\single cell\analysis_output"
STEP = os.path.join(OUT, "step3")
os.makedirs(STEP, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT, ".mplconfig"))

LYO = r"D:\桌面\single cell\DeGene.list.xls"
MARK = os.path.join(OUT, "markers_10groups.csv")
AD = os.path.join(OUT, "step2", "step2_integrated.h5ad")

SETS = {13:"LyO_Progenitor", 4:"LyO_Lymphocyte_precursor", 2:"LyO_Macrophage_like",
        7:"LyO_Hemocytes", 0:"LyO_ECM_producing"}


def load_sets(a):
    d = pd.read_excel(LYO, sheet_name="Top50.DeGene.list")
    import collections
    sym_map = {}
    sym_count = collections.Counter()
    for i, s in enumerate(a.var["sym"].astype(str)):
        sym_map.setdefault(s.lower(), a.var_names[i]); sym_count[s.lower()] += 1
    vn = set(map(str, a.var_names))
    out = {}
    for cl, name in SETS.items():
        g = d[d.Target_Cluster == cl].sort_values("Log2FC", ascending=False).head(50)
        genes, labels = [], []
        for x in g["Gene Name"].astype(str):
            if x in vn:
                genes.append(x); labels.append(x); continue
            if sym_count.get(x.lower(), 0) != 1:
                continue
            v = sym_map.get(x.lower())
            if v is not None and v not in genes:
                genes.append(v); labels.append(x)
        out[name] = {"idx": [list(a.var_names).index(v) for v in genes], "genes": genes, "labels": labels}
    return out


def ucell_scores(a, sets, chunk=300):
    n, N = a.n_obs, a.n_vars
    scores = {k: np.zeros(n, dtype=np.float32) for k in sets}
    for start in range(0, n, chunk):
        end = min(n, start + chunk)
        X = np.asarray(a[start:end].X.todense(), dtype=np.float32)
        order = np.argsort(-X, axis=1, kind="stable")
        ranks = np.empty_like(order)
        rows = np.arange(end - start)[:, None]
        ranks[rows, order] = np.arange(1, N + 1, dtype=np.int32)
        for k, v in sets.items():
            idx = v["idx"]; kk = len(idx)
            if kk == 0: continue
            S = ranks[:, idx].sum(1)
            scores[k][start:end] = 1.0 - (S - kk*(kk+1)/2) / (kk*(N-kk))
    return scores


def main():
    a = ad.read_h5ad(AD)
    sets = load_sets(a)
    for k, v in sets.items():
        print(k, "mapped", len(v["genes"]))
    scores = ucell_scores(a, sets)
    for k, v in scores.items():
        a.obs[f"ucell_{k}"] = v
    a.write_h5ad(os.path.join(STEP, "step3_lyo_scores.h5ad"), compression="gzip")

    order = ["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
    rows = []
    for ct in order:
        m = (a.obs.cell_type == ct).values
        row = {"cell_type":ct,"n":int(m.sum())}
        for k in sets:
            row[f"mean_{k}"] = float(a.obs[f"ucell_{k}"][m].mean())
            row[f"median_{k}"] = float(a.obs[f"ucell_{k}"][m].median())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(STEP,"step3_ucell_by_celltype.csv"), index=False, encoding="utf-8-sig")

    # UMAP panels
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    for ax, k in zip(axes, sets):
        v = a.obs[f"ucell_{k}"].values
        sca = ax.scatter(a.obsm["X_umap"][:,0], a.obsm["X_umap"][:,1], c=v, cmap="viridis", s=2.5,
                         linewidths=0, rasterized=True, vmin=np.quantile(v,0.02), vmax=np.quantile(v,0.98))
        ax.set_title(k, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        plt.colorbar(sca, ax=ax, fraction=0.04)
    axes[-1].axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step3_ucell_umap.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step3_ucell_umap.pdf"),bbox_inches="tight"); plt.close(fig)

    # boxplot of mean scores by cell type
    fig, axes = plt.subplots(1, 2, figsize=(12,4.2))
    for ax, k in zip(axes, ["LyO_Progenitor","LyO_Lymphocyte_precursor"]):
        data = [a.obs.loc[a.obs.cell_type==c, f"ucell_{k}"].values for c in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        ax.set_title(k, fontsize=9); ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step3_ucell_boxplot.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step3_ucell_boxplot.pdf"),bbox_inches="tight"); plt.close(fig)

    # hypergeometric overlap
    mk = pd.read_csv(MARK)
    universe = set(map(str, a.var_names))
    def top_set(grp, n=50):
        g = mk[mk.grp == grp].sort_values("scores", ascending=False).head(n)
        return set(x for x in g["names"].astype(str) if x in universe)
    ph = top_set("PH")
    niche = set().union(*[top_set(g) for g in ["TS","MH2","MyoEpi","M"]])
    rows = []
    for name, v in sets.items():
        A = set(v["genes"])
        for label, B in [("PH_markers", ph), ("niche_markers", niche)] + [(f"{g}_markers", top_set(g)) for g in ["TS","MH2","MyoEpi","M"]]:
            k = len(A & B)
            M = len(universe); n = len(A); N = len(B)
            p = hypergeom.sf(k-1, M, n, N) if k > 0 else 1.0
            exp = n*N/M
            rows.append({"lyo_set":name,"compared_to":label,"overlap":k,"expected":exp,
                         "fold":k/exp if exp>0 else np.nan,"p_hypergeom":p,
                         "overlap_genes":";".join(sorted(A & B))})
    ov = pd.DataFrame(rows).sort_values(["p_hypergeom","overlap"], ascending=[True,False])
    ov.to_csv(os.path.join(STEP,"step3_hypergeom_overlaps.csv"), index=False, encoding="utf-8-sig")

    # heatmap for progenitor set
    fig, ax = plt.subplots(figsize=(7,3.2))
    sub = ov[ov.lyo_set=="LyO_Progenitor"]
    im = ax.imshow(sub[["fold"]].values.T, aspect="auto", cmap="Reds")
    ax.set_xticks(range(len(sub))); ax.set_xticklabels(sub["compared_to"], rotation=45, ha="right")
    ax.set_yticks([0]); ax.set_yticklabels(["fold enrichment"])
    for i, r in enumerate(sub.itertuples()):
        ax.text(i,0,f"{r.overlap} ({r.fold:.1f}x)\np={r.p_hypergeom:.1e}",ha="center",va="center",fontsize=7)
    ax.set_title("LyO Progenitor top50 vs markers", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.03)
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step3_overlap_heatmap.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step3_overlap_heatmap.pdf"),bbox_inches="tight"); plt.close(fig)

    # save gene sets
    pd.DataFrame([{"set":k,"n":len(v["genes"]),"genes":";".join(v["labels"])} for k,v in sets.items()]).to_csv(
        os.path.join(STEP,"step3_genesets.csv"), index=False, encoding="utf-8-sig")
    print(df.to_string(index=False))
    print(ov.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
