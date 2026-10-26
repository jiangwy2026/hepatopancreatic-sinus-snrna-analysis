import os, pandas as pd
from pyscenic.utils import load_motif_annotations
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
for tbl in ("motifs-v9-nr.flybase.tbl", "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"):
    ann = load_motif_annotations(os.path.join(DB, tbl))
    tfs = set(ann.index.get_level_values(0))
    inter = sorted(tfs & set(expr.columns))
    print("%s: annotated TFs %d | present in our %d genes: %d" % (tbl, len(tfs), expr.shape[1], len(inter)))
    print("   sample:", inter[:20])
u = set()
for tbl in ("motifs-v9-nr.flybase.tbl", "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"):
    ann = load_motif_annotations(os.path.join(DB, tbl))
    u |= set(ann.index.get_level_values(0))
print("union annotated TFs:", len(u), "| in our genes:", len(u & set(expr.columns)))
