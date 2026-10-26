import feather, os
import pandas as pd
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
for f in ("dm6_v10_clust.genes_vs_motifs.rankings.feather",):
    t = feather.read_table(os.path.join(DB, f))
    print(f, "rows", t.num_rows, "cols", t.num_columns)
    print("  col names[:6]:", t.column_names[:6])
    print("  last col:", t.column_names[-1])
    ann = pd.read_csv(os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"), sep="\t", low_memory=False)
    db_motifs = set(t.column_names) - {"features", "name"}
    ann_ids = set(ann["#motif_id"].astype(str))
    print("  db motif ids sample:", sorted(db_motifs)[:5])
    print("  ann ids sample:", sorted(ann_ids)[:5])
    print("  overlap db vs ann ids:", len(db_motifs & ann_ids), "of", len(db_motifs))
