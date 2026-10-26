# -*- coding: utf-8 -*-
"""Check which of the manuscript-Table-1 conserved orthologs are detectable in our scRNA-seq,
and in which cell group they are expressed."""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
# genes as named in the manuscript's Table 1 (shrimp gene names)
HS = ["IRF4","neprilysin-1","ferritin","Acsl","cMnSOD3","Tctp","syntenin","STK RIO3","GS","BRI3",
      "Past1","CREG","Arp3","su(r)"]
LO = ["Act5C","Celf2","Cpvl","FAM122A","Mmel1","FCN1","Dsp1","HPGD","Lerp","KCNQ5","LGMN","LYST",
      "mucin2","lymphocyte antigen","Nlrp3","nt5c2","Rgs2","slc16a10","stai","Syne1"]
ALL = HS + LO
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
A.obs["ct"] = pd.Categorical(A.obs.vendor_cluster.astype(str).map(base).replace({"NE": "NEU", "Neu": "NEU"}), categories=ORD)
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym = A.var["sym"].astype(str).values
rows = []
for g in ALL:
    idx = np.where(sym == g)[0]
    if len(idx) == 0:
        rows.append({"gene": g, "found": 0, "best_id": "", "pct_all": 0.0, "top_group": "", "top_pct": 0.0}); continue
    tot = np.asarray(A.layers["counts"][:, idx].sum(0)).ravel()
    j = idx[int(np.argmax(tot))]
    v = np.asarray(A[:, j].X.todense()).ravel()
    pcts = {c: 100 * (v[(A.obs["ct"] == c).values] > 0).mean() for c in ORD}
    top = max(pcts, key=pcts.get)
    rows.append({"gene": g, "found": len(idx), "best_id": A.var_names[j], "pct_all": round(100 * (v > 0).mean(), 2),
                 "top_group": top, "top_pct": round(pcts[top], 1),
                 "pcts": " ".join("%s:%0.0f" % (c, pcts[c]) for c in ORD)})
t = pd.DataFrame(rows)
print(t[["gene", "found", "best_id", "pct_all", "top_group", "top_pct"]].to_string(index=False))
print()
print("----- 各群检出率(%) -----")
for r in rows:
    if r["found"]: print("  %-18s %s" % (r["gene"], r.get("pcts", "")))
t.to_csv(os.path.join(OUT, "manuscript_table1_genes_in_scrna.csv"), index=False)
