# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic import prune as P
from pyscenic.utils import modules_from_adjacencies, load_motif_annotations
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
mtf = set(m.transcription_factor for m in mods)
ann = load_motif_annotations(os.path.join(DB, "motifs-v9-nr.flybase.tbl"))
tfs = set(ann.index.get_level_values(0))
print("module TFs:", len(mtf), "| annotated TFs:", len(tfs))
print("intersection:", len(mtf & tfs))
missing = sorted(mtf - tfs)
print("missing (first 25):", missing[:25])
print("'Awh' in tfs:", "Awh" in tfs, "| 'Awh' in mtf:", "Awh" in mtf)
db = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
sel = [m for m in mods if m.transcription_factor in tfs][:30]
print("testing modules with annotated TF:", len(sel))
df = P.modules2df(db, sel, ann)
print("modules2df ->", df.shape)
if df.shape[0]:
    print(df.sort_values("NES", ascending=False).head(10)[["NES","AUC","Annotation"]].to_string())
