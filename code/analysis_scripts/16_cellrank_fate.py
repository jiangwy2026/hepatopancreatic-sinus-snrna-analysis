import sys, io, os, traceback
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures_traj"
A = sc.read_h5ad(out + r"\adata_hemocyte_lineage.h5ad")
print("cells:", A.shape, "| sub:", A.obs["sub"].value_counts().to_dict())
import cellrank as cr
print("cellrank", cr.__version__)
k = cr.kernels.ConnectivityKernel(A); k.compute_transition_matrix()
est = cr.estimators.GPCCA(k)
est.compute_schur(n_components=10)
est.compute_macrostates(n_states=4, cluster_key="sub")
print("macrostates:", est.macrostates.cat.categories.tolist())
for fn, args in [("set_terminal_states", None), ("predict_terminal_states", None)]:
    pass
try:
    est.predict_terminal_states()
    print("terminal(预测):", list(est.terminal_states.cat.categories))
except Exception as e:
    print("predict_terminal_states 失败:", e)
try:
    est.set_initial_states(["2"], cluster_key="sub")
except TypeError as e:
    print("set_initial_states 签名不同:", e)
    try: est.set_initial_states("2", cluster_key="sub")
    except Exception as e2: print("仍失败:", e2)
except Exception as e:
    print("set_initial_states 失败:", e)
try:
    est.compute_fate_probabilities()
    fp = est.fate_probabilities
    for st in fp.names: A.obs["fate_%s" % st] = np.asarray(fp[st]).ravel()
    cols = ["fate_%s" % s for s in fp.names]
    print("\n=== 命运概率（按原亚群平均）==="); print(A.obs.groupby("sub")[cols].mean().round(3).to_string())
    print("\n=== 命运概率（treated vs blank）==="); print(A.obs.pivot_table(index="sub", columns="group", values=cols, aggfunc="mean").round(3).to_string())
    A.obs[["sub","group","phase","S_score","G2M_score","dpt_pseudotime"]+cols].to_csv(out + r"\fate_probabilities.csv")
    try:
        for lin in fp.names:
            d = est.compute_lineage_drivers(lineages=[lin], cluster_key="sub", return_dataframe=True)
            print("\n=== 命运 %s 的驱动基因 top12 ===" % lin)
            print(d.head(12).round(3).to_string())
    except Exception as e:
        print("lineage drivers 失败:", e)
    sc.settings.figdir = fig
    try: est.plot_fate_probabilities(same_plot=False, show=False)
    except Exception as e: print("plot 失败:", e)
except Exception:
    traceback.print_exc()
A.write(out + r"\adata_hemocyte_lineage.h5ad")
print("\nsaved")
