import os, pyarrow.feather as feather
DB = r"D:\桌面\single cell\analysis_output\scenic_db"
t = feather.read_table(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"))
m = t.column("motifs").to_pylist()
print("motifs sample:", m[:5], "| n:", len(m))
import pandas as pd
ann = pd.read_csv(os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"), sep="\t", low_memory=False)
ids = set(ann["#motif_id"].astype(str))
print("overlap with v10 annotation:", len(set(m) & ids))
ann9 = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
ids9 = set(ann9["#motif_id"].astype(str))
print("overlap with v9 annotation:", len(set(m) & ids9))
