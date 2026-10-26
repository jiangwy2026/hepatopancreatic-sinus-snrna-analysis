import os, pyarrow.feather as feather
DB = r"D:\桌面\single cell\analysis_output\scenic_db"
t = feather.read_table(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"))
print("rows", t.num_rows, "cols", t.num_columns)
print("col names first 8:", t.column_names[:8])
print("col names last 5:", t.column_names[-5:])
import pandas as pd
ann = pd.read_csv(os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"), sep="\t", low_memory=False)
ids = set(ann["#motif_id"].astype(str))
cols = set(t.column_names)
print("annotation ids sample:", sorted(ids)[:4])
print("overlap:", len(cols & ids))
ann9 = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
print("v9 annotation overlap:", len(cols & set(ann9["#motif_id"].astype(str))))
