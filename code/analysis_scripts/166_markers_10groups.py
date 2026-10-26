# -*- coding: utf-8 -*-
"""166_markers_10groups.py - 按最终 10 群重算 marker（含特异度），供热图与单基因图使用"""
import io, sys, os, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
GROUPS = {"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER = ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
MT = {"ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"}
HK = re.compile(r"^(Rp[LS]\d|rpl-|rps-|RpLP|eef|eIF|Act\d|Act5C|Tub|Hsp\d|hsc\d|HSC70|PRSS)", re.I)
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.obs["ct"] = pd.Categorical(A.obs["vendor_cluster"].astype(str).map(GROUPS), categories=ORDER)
A.var["is_mt"] = A.var["sym"].astype(str).isin(MT)
sc.tl.rank_genes_groups(A, "ct", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(A, group=g).assign(grp=g) for g in ORDER], ignore_index=True)
mk["sym"] = mk["names"].map(dict(zip(A.var_names.astype(str), A.var["sym"].astype(str)))).fillna(mk["names"]).astype(str)
mk["is_mt"] = mk["names"].map(A.var["is_mt"].to_dict()).fillna(False)
mk["hk"] = mk["sym"].str.match(HK).fillna(False)
mk["dPct"] = (mk.pct_nz_group - mk.pct_nz_reference) * 100
mk.to_csv(os.path.join(OUT, "markers_10groups.csv"), index=False)
sig = mk[(~mk.is_mt) & (~mk.hk) & (mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)].copy()
print("每组显著 marker:", sig.groupby("grp").size().to_dict())
print("\n=== 每组 top 8（按 score，去重）===")
rows = []; seen = set()
for g in ORDER:
    d = sig[sig.grp == g].sort_values(["scores"], ascending=False)
    show = []
    for r in d.itertuples():
        if r.sym in seen: continue
        seen.add(r.sym); show.append((r.sym, r.names, round(float(r.logfoldchanges), 2), round(float(r.dPct), 1)))
        if len(show) >= 8: break
    print("  %-8s %s" % (g, ", ".join(x[0] for x in show)))
    for rank, (s, n, lfc, dp) in enumerate(show, 1):
        rows.append({"grp": g, "rank": rank, "sym": s, "gid": n, "log2FC": lfc, "dPct": dp})
T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, "top_markers_10groups.csv"), index=False)
print("\nsaved markers_10groups.csv / top_markers_10groups.csv")
