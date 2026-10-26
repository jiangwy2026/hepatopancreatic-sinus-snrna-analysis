# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd, pyarrow.feather as feather
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic import prune as P
from pyscenic.utils import load_motif_annotations, modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")

for f in ("dm6_v10_clust.regions_vs_motifs.rankings.feather",):
    t = feather.read_table(os.path.join(DB, f))
    last = t.column_names[-1]
    m = set(t.column(last).to_pylist())
    for tbl in ("motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl", "motifs-v9-nr.flybase.tbl"):
        ann = pd.read_csv(os.path.join(DB, tbl), sep="\t", low_memory=False)
        ids = set(ann["#motif_id"].astype(str))
        print("%s (rows %d, last col %s): overlap with %s = %d / %d" % (
            f, t.num_rows, last, tbl, len(m & ids), len(m)))

expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
db = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
for tbl in ("motifs-v9-nr.flybase.tbl", "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"):
    ann = load_motif_annotations(os.path.join(DB, tbl))
    df = P.modules2df(db, mods[:5], ann)
    print("--- %s -> modules2df shape %s" % (tbl, df.shape))
    if df.shape[0]:
        print(df.sort_values("NES", ascending=False).head(6)[["NES","AUC","Annotation"]].to_string())
        break
