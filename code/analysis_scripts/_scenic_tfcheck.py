import os, numpy as np, pandas as pd
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic.utils import modules_from_adjacencies, load_motif_annotations
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
mtf = [m.transcription_factor for m in mods]
print("module TFs (n=%d):" % len(mtf), mtf[:15])
for tbl in ("motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl", "motifs-v9-nr.flybase.tbl"):
    ann = load_motif_annotations(os.path.join(DB, tbl))
    tfs = set(ann.index.get_level_values(0))
    awh = ann.loc[ann.index.get_level_values(0).isin(mtf)]
    print("=== %s: %d rows | unique TFs %d | module TFs with any annotation: %d/%d" % (
        tbl, len(ann), len(tfs), len(set(mtf) & tfs), len(mtf)))
