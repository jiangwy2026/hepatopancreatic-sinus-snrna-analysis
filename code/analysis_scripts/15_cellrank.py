import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_hemocyte_dpt.h5ad")
A = A[A.obs["sub"].isin(["0","1","2"])].copy()
print("细胞数:", A.shape[0])
sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
import cellrank as cr
print("cellrank", cr.__version__)
k = cr.kernels.ConnectivityKernel(A)
k.compute_transition_matrix()
est = cr.estimators.GPCCA(k)
print("estimator methods:", [m for m in dir(est) if not m.startswith("_")][:40])
try:
    est.compute_schur(n_components=10)
except Exception as e:
    print("schur fail:", e); est.compute_eigendecomposition(k=10)
est.compute_macrostates(n_states=5, cluster_key="sub")
print("\nmacrostates:")
print(est.macrostates.cat.categories.tolist())
try:
    est.set_initial_states(states=["2"], cluster_key="sub")
    print("initial set ->", est.initial_states)
except Exception as e:
    print("set_initial_states fail:", e)
try:
    est.predict_terminal_states()
    print("terminal ->", est.terminal_states)
except Exception as e:
    print("predict_terminal fail:", e)
try:
    est.compute_fate_probabilities()
    fp = est.fate_probabilities
    print("\n=== 命运概率（按原亚群平均）===")
    for st in fp.names:
        A.obs["fate_" + str(st)] = np.asarray(fp[st]).ravel()
    cols = ["fate_" + str(s) for s in fp.names]
    print(A.obs.groupby("sub")[cols].mean().round(3).to_string())
    print("\n=== 命运概率：treated vs blank（按原亚群）===")
    print(A.obs.pivot_table(index="sub", columns="group", values=cols, aggfunc="mean").round(3).to_string())
    A.obs[["sub","group","phase","G2M_score","S_score"]+cols].to_csv(out + r"\fate_probabilities.csv")
except Exception as e:
    import traceback; traceback.print_exc()
A.write(out + r"\adata_hemocyte_lineage.h5ad")
print("\nsaved")
