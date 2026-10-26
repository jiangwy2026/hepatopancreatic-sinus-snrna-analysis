import sys, io, os, time
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
base = {"0":"MH1","1":"GH","2":"PH","3":"MH2"}
A.obs["ct"] = A.obs.vendor_cluster.astype(str).map(base)
H = A[A.obs["ct"].notna()].copy()
H.obs["ct"] = pd.Categorical(H.obs["ct"], categories=["PH","GH","MH1","MH2"])
print("血细胞（含 MH2）:", H.shape, dict(H.obs["ct"].value_counts()))
# 用校正后的 PCA 空间（重建 Harmony 坐标：从 hvg 对象取）
Ah = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")
idx = [list(Ah.obs_names).index(x) for x in H.obs_names]
H.obsm["X_pca_harmony"] = Ah.obsm["X_pca_harmony"][idx]
sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
# 扩散映射
sc.tl.diffmap(H, n_comps=15)
dc = H.obsm["X_diffmap"]
print("DC1 范围 %.3f ~ %.3f | DC2 %.3f ~ %.3f" % (dc[:,1].min(), dc[:,1].max(), dc[:,2].min(), dc[:,2].max()))
# root = PH 群中 PH 定位最强的细胞（用 PH marker）
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)): sym2i.setdefault(s, i)
phmk = [g for g in ["LOC113823945","ccnd2","Igfbp7","Grik2"] if g in sym2i]
score = np.asarray(H[:, [sym2i[g] for g in phmk]].X.mean(1)).ravel()
isph = (H.obs["ct"] == "PH").values
iroot = int(np.where(isph)[0][np.argmax(score[isph])])
H.uns["iroot"] = iroot
sc.tl.dpt(H)
pt = H.obs["dpt_pseudotime"].values.astype(float)
print("root 细胞:", H.obs_names[iroot], "| 拟时序 有限=%d/%d | 范围 %.3f~%.3f" % (np.isfinite(pt).sum(), len(pt), np.nanmin(pt), np.nanmax(pt)))
import pickle; pickle.dump({"iroot": iroot}, open(out + r"\_fig5_root.pkl", "wb"))
H.write(out + r"\adata_hemocyte_fig5.h5ad")
print("saved adata_hemocyte_fig5.h5ad")
