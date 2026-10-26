# -*- coding: utf-8 -*-
"""230_step1_niche_programs.py

Step 1a: define candidate niche cells in the hepatic-sinusoid dataset and
quantify niche-support programmes.

Input : analysis_output/adata_soupx.h5ad      (9658 cells, raw counts in layers['counts'])
        analysis_output/gene_annotation.csv   (NCBI product names per GeneID)
Output: analysis_output/step1/...
"""
import os, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\桌面\single cell\analysis_output"
STEP = os.path.join(OUT, "step1")
os.makedirs(STEP, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT, ".mplconfig"))

GROUP = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
         "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
NICHE_COMPARTMENTS = ["TS", "MH2", "MyoEpi", "M"]
HEMOCYTE = ["GH", "PH"]
HEPATOPANCREAS = ["F", "R", "B"]

# Gene panels defined by NCBI GeneID (product-name checked), not by the
# vendor symbol column, because many symbols are mis-assigned in this annotation.
PANEL = {
 "Endothelial": [("113816763","TIE1"),("113814550","Tie1-like"),("113803203","Tie1-like/V-RYK"),
                 ("113803634","VEGFR2-like"),("113806166","VEGFR1-like"),("113811468","VEGFR3-like"),
                 ("113821981","PTPR-alpha-like")],
 "Stromal_ECM": [("113824160","COMP/THBS4-like"),("113830308","COMP-like"),("113803248","THBS4-B-like"),
                 ("113808921","Laminin beta-1"),("113814443","Laminin gamma-1"),("113823054","Laminin alpha-1"),
                 ("113804315","MMP2-like"),("113808646","MMP14-like"),("113803924","Collagen XVIII")],
 "Macrophage": [("113822909","MRC1-like"),("113805228","Mannose receptor"),("113824479","PDGFR-alpha/CSF1R-like"),
                ("113810024","MH2 chitin-binding"),("113813946","MH2 marker"),("113811918","MH2 marker"),
                ("113821858","Jagged-1b")],
 "Hematopoietic_support": [("113822052","Astakine"),("113805969","PPAF2"),("113809349","PPAF1"),
                ("113822272","PPAF3"),("113822289","PPO-activating enzyme"),("113824171","Hemocytin"),
                ("113821858","Jagged-1b"),("113815864","Delta-like 1")],
 "Hemocyte_progenitor": [("113805969","PPAF2"),("113809687","A2M"),("113824171","Hemocytin"),
                ("113812814","PPAF2-like")],
 "Hepatopancreas": [("113810052","Trypsin-like"),("113822654","Chymotrypsin-like"),("113818984","ACE-like"),
                ("113805758","Trypsin-3-like"),("113817089","Jon99Cii-like")],
}
# receptor / target panel used later by the communication script
RECEPTORS = {
 "Notch3-like":"113812198", "Notch-like":"113802681", "Notch2-like":"113812690", "Notch4-like":"113812693",
 "Integrin alpha-V":"113802877", "Integrin alpha-5":"113807444", "Integrin beta-PS":"113805715",
 "TIE1":"113816763", "Tie1-like":"113814550", "IGF1R":"113821898",
}


def main():
    a = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    a.X = a.layers["counts"].copy()
    ann = pd.read_csv(os.path.join(OUT, "gene_annotation.csv"), dtype={"GeneID": str})[["GeneID","name"]]
    old_names = [str(x) for x in a.var_names]
    v = a.var.copy()
    v["gid"] = [str(x).replace("ncbi_", "") for x in old_names]
    v = v.merge(ann, left_on="gid", right_on="GeneID", how="left", suffixes=("","_ann"))
    v.index = old_names
    a.var = v
    gid2idx = {g: i for i, g in enumerate(a.var["gid"].astype(str))}
    gid2var = {g: a.var_names[i] for g, i in gid2idx.items()}

    a.obs["cell_type"] = a.obs["vendor_cluster"].astype(str).map(GROUP)
    a.obs["compartment"] = a.obs["cell_type"].map(lambda x: "niche" if x in NICHE_COMPARTMENTS else
                                                  ("hemocyte" if x in HEMOCYTE else
                                                   ("hepatopancreas" if x in HEPATOPANCREAS else "other")))

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers["lognorm"] = a.X.copy()

    # module scores
    score_rows = []
    for mod, pairs in PANEL.items():
        genes = [gid2var[g] for g, _ in pairs if g in gid2var]
        if len(genes) < 2:
            continue
        key = f"score_{mod}"
        sc.tl.score_genes(a, genes, score_name=key, ctrl_size=50, random_state=0)
        score_rows.append({"module": mod, "n_genes": len(genes), "genes": ";".join(genes)})
    pd.DataFrame(score_rows).to_csv(os.path.join(STEP, "step1_scored_modules.csv"), index=False, encoding="utf-8-sig")

    # niche subtype per cell (only within niche compartments), z-scored so a
    # globally high but non-specific module cannot dominate the argmax
    mod_cols = {"Endothelial_like":"z_Endothelial", "Stromal_like":"z_Stromal_ECM",
                "Macrophage_like":"z_Macrophage"}
    for mod in ["Endothelial","Stromal_ECM","Macrophage","Hematopoietic_support"]:
        v = a.obs[f"score_{mod}"].astype(float)
        a.obs[f"z_{mod}"] = (v - v.mean()) / (v.std() + 1e-9)
    sub = pd.Series("non-niche", index=a.obs_names, dtype=object)
    niche = a.obs["compartment"].eq("niche").values
    S = a.obs[list(mod_cols.values())].values
    best = np.argmax(S, axis=1)
    names = np.array(list(mod_cols.keys()))
    sub.iloc[np.where(niche)[0]] = names[best[niche]]
    a.obs["niche_subtype"] = sub.values

    # embedding (batch-aware HVG)
    b = a.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, batch_key="sample")
    b = b[:, b.var["highly_variable"].values].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=50, svd_solver="arpack", random_state=0)
    sc.pp.neighbors(b, n_neighbors=15, n_pcs=40, random_state=0)
    sc.tl.umap(b, random_state=0)
    a.obsm["X_umap"] = b.obsm["X_umap"]
    a.write_h5ad(os.path.join(STEP, "step1_niche.h5ad"), compression="gzip")

    cols = ["sample","group","cell_type","compartment","niche_subtype","vendor_cluster"] + list(mod_cols.values())
    a.obs[cols].to_csv(os.path.join(STEP, "step1_cell_scores.csv"), encoding="utf-8-sig")

    # composition
    comp = pd.crosstab(a.obs["cell_type"], a.obs["niche_subtype"], normalize="index")
    comp.to_csv(os.path.join(STEP, "step1_niche_subtype_composition.csv"), encoding="utf-8-sig")

    # ---- Figures (English labels, publication-ish) ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, key, title in zip(axes, ["cell_type","compartment","niche_subtype"],
                              ["Cell type","Compartment","Niche subtype"]):
        cats = a.obs[key].astype(str).unique().tolist()
        cat = a.obs[key].astype("category").cat.set_categories(cats)
        ax.scatter(a.obsm["X_umap"][:,0], a.obsm["X_umap"][:,1], c=cat.cat.codes,
                   cmap="tab20", s=3, linewidths=0, rasterized=True)
        for i, c in enumerate(cats):
            m = cat == c
            if m.sum() == 0: continue
            ax.text(np.median(a.obsm["X_umap"][m.values,0]), np.median(a.obsm["X_umap"][m.values,1]),
                    c, fontsize=7, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.75))
        ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(STEP, "step1_umap.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(STEP, "step1_umap.pdf"), bbox_inches="tight")
    plt.close(fig)

    # feature plots of module scores
    mods = list(mod_cols.values())
    fig, axes = plt.subplots(1, len(mods), figsize=(4.2*len(mods), 4.0))
    for ax, m in zip(np.atleast_1d(axes), mods):
        v = a.obs[m].values
        sca = ax.scatter(a.obsm["X_umap"][:,0], a.obsm["X_umap"][:,1], c=v, cmap="Reds", s=3,
                         linewidths=0, rasterized=True, vmin=np.quantile(v,0.01), vmax=np.quantile(v,0.99))
        ax.set_title(m.replace("score_","").replace("_"," "), fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        plt.colorbar(sca, ax=ax, fraction=0.04)
    fig.tight_layout()
    fig.savefig(os.path.join(STEP, "step1_module_scores_umap.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(STEP, "step1_module_scores_umap.pdf"), bbox_inches="tight")
    plt.close(fig)

    # dotplot of niche factors across cell types (average lognorm expression + pct)
    factor_ids = []
    for mod in ["Endothelial","Stromal_ECM","Macrophage","Hematopoietic_support"]:
        for g, lab in PANEL[mod]:
            if g in gid2var and (g, lab) not in factor_ids:
                factor_ids.append((g, lab))
    genes = [gid2var[g] for g, _ in factor_ids if g in gid2var]
    labels = [lab for g, lab in factor_ids if g in gid2var]
    seen = {}
    uniq = []
    for lab in labels:
        seen[lab] = seen.get(lab, 0) + 1
        uniq.append(lab if seen[lab] == 1 else f"{lab} ({seen[lab]})")
    labels = uniq
    d = a[:, genes].copy()
    d.X = d.layers["lognorm"].copy()
    d.var["label"] = labels
    d.var_names = labels
    sc.pl.dotplot(d, labels, groupby="cell_type", standard_scale="var", show=False,
                  cmap="Reds", figsize=(max(6, 0.34*len(labels)), 3.2))
    plt.savefig(os.path.join(STEP, "step1_niche_factor_dotplot.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(STEP, "step1_niche_factor_dotplot.pdf"), bbox_inches="tight")
    plt.close("all")

    # group-level mean expression table for the same panel
    rows = []
    for g, lab in factor_ids:
        if g not in gid2var: continue
        idx = gid2var[g]
        x = np.asarray(a[:, idx].layers["lognorm"].toarray()).ravel()
        row = {"gene": lab, "gid": g}
        for ct in ["TS","GH","PH","MH2","PM-like","F","R","B","MyoEpi","M"]:
            m = a.obs["cell_type"].eq(ct).values
            row[ct] = float(x[m].mean()) if m.sum() else np.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(STEP, "step1_niche_factor_group_means.csv"), index=False, encoding="utf-8-sig")

    # summary json
    summary = {
        "n_cells": int(a.n_obs),
        "cell_type_counts": {k: int(v) for k, v in a.obs["cell_type"].value_counts().items()},
        "niche_subtype_counts": {k: int(v) for k, v in a.obs["niche_subtype"].value_counts().items()},
        "niche_compartments": NICHE_COMPARTMENTS,
        "scored_modules": {m: len([1 for g,_ in p if g in gid2var]) for m,p in PANEL.items()},
    }
    with open(os.path.join(STEP, "step1_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
