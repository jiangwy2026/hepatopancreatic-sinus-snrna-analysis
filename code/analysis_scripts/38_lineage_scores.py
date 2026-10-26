import sys, io, os, re
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 300)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
HEM = ["A2M","PPAF1","PPAF2","Pxt","SSPO","TGM1","LOC113823945","LOC113801825","LOC113816267","PEN-3","LOC113800363","LOC113813598","FCN2","SERPINB1"]
HP  = ["Prss1","Prss2","Try3","Jon99Cii","Cpa2","LCP2","Chia","ACE","LIPF","Mal-A2","MALRD1","Oacyl","HSPG2","LOC113810052"]
CC  = ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","CCNB3","AURKA","TOP2A","TPX2","KIF20A","KIF23","RACGAP1","INCENP","ASPM","Pka-C1","Patronin"]
sets = {"血细胞身份":HEM, "肝胰腺身份":HP, "细胞周期":CC}
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
id2sym = {num(g): s for g, s in zip(feats["gid"], feats["sym"])}
def s2i_build(A):
    if "sym" in A.var.columns:
        d = {}
        for i, s in zip(A.var_names, A.var["sym"]):
            if s is not None and pd.notna(s): d.setdefault(str(s), i)
        return d
    return {id2sym[k]: i for i, k in enumerate(A.var_names) if k in id2sym}
def score_ds(A, label, groupcol):
    A = A.copy(); A.X = A.layers["counts"].copy() if "counts" in A.layers else A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    s2i = s2i_build(A)
    labs = [c for c in set(A.obs[groupcol].astype(str)) if str(c).isdigit()]
    rows = []
    for c in sorted(labs, key=int):
        m = (A.obs[groupcol].astype(str) == c).values
        if m.sum() < 20: continue
        r = {"dataset": label, "cluster": c, "n": int(m.sum())}
        for nm, gs in sets.items():
            ids = [s2i[g] for g in gs if g in s2i]
            r[nm] = float(np.asarray(A[m][:, ids].X.mean(0)).ravel().mean()) if len(ids) >= 3 else np.nan
            r[nm + "_ngene"] = len(ids)
        rows.append(r)
    return pd.DataFrame(rows)
M = sc.read_h5ad(out + r"\adata_soupx.h5ad"); M.var["sym"] = feats.set_index("gid").reindex(M.var_names)["sym"].values
d1 = score_ds(M, "我的数据", "vendor_cluster")
E = sc.read_h5ad(out + r"\elife_adata.h5ad"); d2 = score_ds(E, "eLife血细胞", "elife_cluster")
G = sc.read_h5ad(out + r"\geo_GSE306674_processed.h5ad")
G.var["sym"] = [id2sym.get(str(x), None) for x in G.var_names]
d3 = score_ds(G, "GEO肝胰腺", "cl")
for d, t in ((d1,"我的数据（10 群）"), (d2,"eLife 血细胞（13 群）"), (d3,"GEO 肝胰腺（17 群）")):
    print("\n=== %s ===" % t)
    cols = ["cluster","n","血细胞身份","肝胰腺身份","细胞周期"]
    print(d[cols].round(3).to_string(index=False))
    print("  基因数: 血细胞 %d | 肝胰腺 %d | 细胞周期 %d" % (d["血细胞身份_ngene"].iloc[0], d["肝胰腺身份_ngene"].iloc[0], d["细胞周期_ngene"].iloc[0]))
pd.concat([d1,d2,d3]).to_csv(out + r"\lineage_scores_all.csv", index=False)
print("\nsaved lineage_scores_all.csv")
