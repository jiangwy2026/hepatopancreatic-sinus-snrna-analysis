import sys, io, os, re
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 300)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; REF = out + r"\elife_ref"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
id2sym = {num(g): s for g, s in zip(feats["gid"], feats["sym"])}
sets = {}
for f, lab in (("elife-80127-supp3-v2.csv","PH1"),("elife-80127-supp6-v2.csv","PH2"),
               ("elife-80127-supp5-v2.csv","MH2"),("elife-80127-supp8-v2.csv","MH1"),
               ("elife-80127-supp4-v2.csv","GH2"),("elife-80127-supp7-v2.csv","GH1")):
    d = pd.read_csv(os.path.join(REF,f), encoding="utf-8-sig"); d = d[d.Qvalue<0.05]
    sets[lab] = d["Gene ID"].astype(str).map(num).unique().tolist()
CYCPAT = re.compile(r"(?i)^(ccn|cdk|cdc|mcm|orc|pola|prim|fen1|lig1|tyms|dhrf|gins|pcna|aur|plk|bub|mad|top2|tpx2|kif|cenp|smc|anln|racgap|ect2|cks|tubb|tuba|espl|ndc80|nuf2|incenp|prc1|aspm|mki67|ube2|h1-|h2a|h2b|h3|h4|rpl|rps|eef|act|tub)")
sets = {k: [x for x in v if not CYCPAT.match(str(id2sym.get(x, x)))] for k, v in sets.items()}
for k in sets: print("%-4s 剔除周期后 %4d 个基因" % (k, len(sets[k])))
def tab(A, label, groupcol):
    A = A.copy(); A.X = A.layers["counts"].copy() if "counts" in A.layers else A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    s2i = {}
    if "sym" in A.var.columns:
        for i, s in zip(A.var_names, A.var["sym"]):
            if s is not None and pd.notna(s): s2i.setdefault(str(s), i)
    else:
        s2i = {id2sym[k2]: i for i, k2 in enumerate(A.var_names) if k2 in id2sym}
    labs = [c for c in set(A.obs[groupcol].astype(str)) if str(c).isdigit()]
    rows = {}
    for c in sorted(labs, key=int):
        m = (A.obs[groupcol].astype(str) == c).values
        if m.sum() < 20: continue
        for k, ids in sets.items():
            idx = sorted({s2i[id2sym.get(x,x)] for x in ids if id2sym.get(x,x) in s2i})
            rows.setdefault(c, {})[k] = float(np.asarray(A[m][:, idx].X.mean(0)).ravel().mean())
    T = pd.DataFrame(rows).T
    Z = T.sub(T.mean(0), axis=1).div(T.std(0) + 1e-9, axis=1)
    return Z
M  = sc.read_h5ad(out + r"\adata_soupx.h5ad"); M.var["sym"] = feats.set_index("gid").reindex(M.var_names)["sym"].values
Mu = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str)
Mu = Mu[Mu.obs_names.isin(set(cl["Cells"]))].copy(); Mu.obs["vendor_cluster"] = cl.set_index("Cells").reindex(Mu.obs_names)["Cluster"].astype(str).values
Mu.var["sym"] = feats.set_index("gid").reindex(Mu.var_names)["sym"].values
E  = sc.read_h5ad(out + r"\elife_adata.h5ad")
G  = sc.read_h5ad(out + r"\geo_GSE306674_processed.h5ad"); G.var["sym"] = [id2sym.get(str(x), None) for x in G.var_names]
for Z, name in ((tab(M,"我的(校正后)","vendor_cluster"), "我的数据·SoupX校正后"),
                (tab(Mu,"我的(未校正)","vendor_cluster"), "我的数据·未校正"),
                (tab(E,"eLife","elife_cluster"), "eLife 血细胞(正对照)"),
                (tab(G,"GEO","cl"), "GEO 肝胰腺")):
    print("\n########## %s（行内 z-score，剔除细胞周期基因后）" % name)
    print(Z.round(2).to_string())
print("\n=== 各数据集中 PH1+PH2 综合得分最高的群 ===")
for Z, name in ((tab(M,"","vendor_cluster"), "我的·校正后"), (tab(Mu,"","vendor_cluster"), "我的·未校正"),
                (tab(E,"","elife_cluster"), "eLife(正对照)"), (tab(G,"","cl"), "GEO 肝胰腺")):
    Z["PH综合"] = Z["PH1"] + Z["PH2"]
    o = Z["PH综合"].sort_values(ascending=False)
    print("  %-16s 最高: %s (%.2f) | 次高: %s (%.2f)" % (name, o.index[0], o.iloc[0], o.index[1], o.iloc[1]))
