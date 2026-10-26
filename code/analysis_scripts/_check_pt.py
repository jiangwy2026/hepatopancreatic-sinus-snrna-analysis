import anndata as ad, numpy as np, os
OUT = r"D:\桌面\single cell\analysis_output"
for f in ("adata_hemocyte_traj.h5ad", "adata_hemocyte_fig4.h5ad", "adata_hemocyte_fig5c.h5ad", "adata_hemocyte_fig5.h5ad"):
    p = os.path.join(OUT, f)
    if not os.path.exists(p):
        print(f, "missing"); continue
    A = ad.read_h5ad(p, backed="r")
    print("==", f, A.shape)
    print("   obs:", list(A.obs.columns))
    print("   obsm:", list(A.obsm.keys()), "| uns:", list(A.uns.keys())[:12])
    for c in A.obs.columns:
        if "time" in c.lower() or "dpt" in c.lower() or "pseudo" in c.lower():
            v = np.asarray(A.obs[c], dtype=float)
            print("   %s: min %.3f p25 %.3f median %.3f p75 %.3f max %.3f" % (
                c, np.nanmin(v), np.nanpercentile(v, 25), np.nanmedian(v), np.nanpercentile(v, 75), np.nanmax(v)))
    if "diffmap_evals" in A.uns:
        print("   diffmap_evals:", np.round(np.asarray(A.uns["diffmap_evals"]), 4)[:10])
