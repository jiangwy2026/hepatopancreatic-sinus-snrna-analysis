import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
def mix(H):
    K = H.obsp["connectivities"].tocoo(); d = H.obs["dataset"].values
    return 100*(1 - (d[K.row] == d[K.col]).mean())
run = []
for label, keep in (("全部细胞", ["0","1","2","3","4","5","6","7","8","9"]),
                    ("仅血细胞群 C0/C1/C2/C3", ["0","1","2","3"]),
                    ("仅 C1+C2（颗粒+增殖）", ["1","2"])):
    A = sc.read_h5ad(out + r"\adata_raw.h5ad")
    cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
    A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
    A.obs["my_cluster"] = cl.set_index("Cells").reindex(A.obs_names)["Cluster"].astype(str).values
    A = A[A.obs.my_cluster.isin(keep)].copy(); A.obs["dataset"]="Mine"
    E = sc.read_h5ad(out + r"\elife_adata.h5ad"); E.obs["dataset"]="eLife"
    shared = sorted(set(A.var_names) & set(E.var_names))
    C = ad.concat([A[:,shared].copy(), E[:,shared].copy()], join="inner", index_unique=None)
    C.obs["dataset"]=C.obs["dataset"].astype(str); C.layers["counts"]=C.X.copy()
    C.X = C.layers["counts"].copy(); sc.pp.normalize_total(C, target_sum=1e4); sc.pp.log1p(C)
    sc.pp.highly_variable_genes(C, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    B = C[:, C.var.highly_variable].copy(); B.X = B.layers["counts"].copy()
    sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B); sc.pp.scale(B, max_value=10)
    sc.tl.pca(B, n_comps=50, svd_solver="arpack")
    m0 = None
    sc.external.pp.harmony_integrate(B, "dataset", basis="X_pca", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(B, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
    m1 = mix(B)
    B2 = B.copy(); sc.pp.neighbors(B2, use_rep="X_pca", n_neighbors=15, n_pcs=30)
    m0 = mix(B2)
    print("%-28s 细胞 %6d (我 %5d / eLife %5d) | 未整合 %.1f%% | Harmony %.1f%%" %
          (label, B.n_obs, (B.obs.dataset=="Mine").sum(), (B.obs.dataset=="eLife").sum(), m0, m1), flush=True)
    run.append(dict(subset=label, n_mine=int((B.obs.dataset=="Mine").sum()), n_elife=int((B.obs.dataset=="eLife").sum()),
                    mix_raw=m0, mix_harmony=m1))
    if label.startswith("仅血细胞"):
        sc.tl.umap(B); um=B.obsm["X_umap"]; d=B.obs.dataset.values
        fig, ax = plt.subplots(figsize=(6,5.2))
        for name, col in (("eLife","#4C72B0"),("Mine","#C44E52")):
            m = d==name; ax.scatter(um[m,0],um[m,1],s=0.7,c=col,linewidths=0,rasterized=True,label="%s (n=%d)"%(name,m.sum()))
        ax.legend(markerscale=8,fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title("Hemocyte compartment only: Mine vs eLife", fontsize=10)
        plt.tight_layout(); plt.savefig(out + r"\figures_joint\joint_umap_hemocytes.png", dpi=150); plt.close()
        B.write(out + r"\joint_hemocytes_only.h5ad")
pd.DataFrame(run).to_csv(out + r"\joint_mixing_summary.csv", index=False)
print("\nsaved joint_mixing_summary.csv")
