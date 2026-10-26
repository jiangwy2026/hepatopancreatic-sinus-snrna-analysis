import sys, io, os, re
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 300)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
id2sym = {num(g): s for g, s in zip(feats["gid"], feats["sym"])}
# 只保留“谱系独占”的基因（互斥）
HEM_ONLY = ["TGM1","LOC113823945","A2M","PPAF1","PPAF2","SSPO","LOC113801825","LOC113816267","PEN-3","LOC113800363","Pxt"]
HP_ONLY  = ["Prss1","Prss2","Try3","Jon99Cii","Cpa2","LCP2","Chia","ACE","LIPF","Mal-A2","MALRD1","LOC113810052","Oacyl"]
CCSET    = ["PCNA","MCM6","CDC6","ORC1","POLA2","PRIM1","FEN1","TYMS","CCNB3","AURKA","TOP2A","TPX2","KIF20A","KIF23","RACGAP1","ASPM","Pka-C1","Patronin"]
def got(A):
    if "sym" in A.var.columns:
        d = {}
        for i, s in zip(A.var_names, A.var["sym"]):
            if s is not None and pd.notna(s): d.setdefault(str(s), i)
        return d
    return {id2sym[k]: i for i, k in enumerate(A.var_names) if k in id2sym}
def run(A, label, groupcol, keep=None):
    A = A.copy(); A.X = A.layers["counts"].copy() if "counts" in A.layers else A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    s2i = got(A)
    def sc_(gs): 
        idx = sorted({s2i[g] for g in gs if g in s2i})
        return idx, len(idx)
    ih, nH = sc_(HEM_ONLY); ip, nP = sc_(HP_ONLY); ic, nC = sc_(CCSET)
    labs = [c for c in set(A.obs[groupcol].astype(str)) if str(c).isdigit()]
    rows = []
    for c in sorted(labs, key=int):
        m = (A.obs[groupcol].astype(str) == c).values
        if m.sum() < 20: continue
        h = float(np.asarray(A[m][:, ih].X.mean(0)).ravel().mean())
        p = float(np.asarray(A[m][:, ip].X.mean(0)).ravel().mean())
        cc = float(np.asarray(A[m][:, ic].X.mean(0)).ravel().mean())
        rows.append({"dataset": label, "cluster": c, "n": int(m.sum()), "血细胞": round(h,3), "肝胰腺": round(p,3),
                     "细胞周期": round(cc,3), "血/肝比值": round(h/(p+0.01), 1)})
    print("\n=== %s (可用基因: 血细胞 %d | 肝胰腺 %d | 周期 %d) ===" % (label, nH, nP, nC))
    d = pd.DataFrame(rows); print(d.to_string(index=False)); return d
M  = sc.read_h5ad(out + r"\adata_soupx.h5ad"); M.var["sym"] = feats.set_index("gid").reindex(M.var_names)["sym"].values
Mu = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
Mu = Mu[Mu.obs_names.isin(set(cl["Cells"]))].copy(); Mu.obs["vendor_cluster"] = cl.set_index("Cells").reindex(Mu.obs_names)["Cluster"].astype(str).values
Mu.var["sym"] = feats.set_index("gid").reindex(Mu.var_names)["sym"].values
E  = sc.read_h5ad(out + r"\elife_adata.h5ad")
G  = sc.read_h5ad(out + r"\geo_GSE306674_processed.h5ad"); G.var["sym"] = [id2sym.get(str(x), None) for x in G.var_names]
a = run(M, "我的数据（SoupX 校正后）", "vendor_cluster")
b = run(Mu, "我的数据（未校正）", "vendor_cluster")
c = run(E, "eLife 血淋巴（正对照）", "elife_cluster")
d = run(G, "GEO 肝胰腺", "cl")
pd.concat([a,b,c,d]).to_csv(out + r"\lineage_ratio_table.csv", index=False)
print("\nsaved lineage_ratio_table.csv")
