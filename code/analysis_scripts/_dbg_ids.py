import io, sys, os
import numpy as np, pandas as pd, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
print("raw.var_names[:3]:", list(raw.var_names[:3]))
print("raw.var 列:", list(raw.var.columns))
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"))
print("four['gid'][:3]:", four["gid"].tolist()[:3], "| 类型:", type(four["gid"].iloc[0]))
gp = {g: i for i, g in enumerate(raw.var_names.astype(str))}
hit = sum(1 for g in four["gid"].astype(str) if g in gp)
print("four gid 命中 raw.var_names: %d/%d" % (hit, len(four)))
nine = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected_fixed.csv"))
print("nine 列:", list(nine.columns)[:8])
print("nine['names'][:3]:", nine["names"].astype(str).tolist()[:3])
print("nine names 命中: %d/%d" % (sum(1 for g in nine["names"].astype(str) if g in gp), len(nine)))
