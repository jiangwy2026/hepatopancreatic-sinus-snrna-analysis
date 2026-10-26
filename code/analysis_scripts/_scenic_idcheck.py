import os, pandas as pd
DB = r"D:\桌面\single cell\analysis_output\scenic_db"
tests = ["metacluster_4.3", "stark__YGATAC", "nitta__prd_TACGGA20NCG_KA_TAATCRATTA_m1_c5", "bergman__Adf1"]
for tbl in ("motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl", "motifs-v9-nr.flybase.tbl"):
    ann = pd.read_csv(os.path.join(DB, tbl), sep="\t", low_memory=False)
    ids = set(ann["#motif_id"].astype(str))
    pref = pd.Series([i.split("__")[0] for i in ids]).value_counts()
    print("===", tbl, "rows", len(ann))
    for t in tests:
        print("   ", t, "->", t in ids)
    print("   prefix counts:", pref.head(10).to_dict())
