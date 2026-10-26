import sys, io, os, re
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 300)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; REF = out + r"\elife_ref"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
id2sym = {num(g): s for g, s in zip(feats["gid"], feats["sym"])}
# eLife PH1/PH2 marker
sets = {}
for f, lab in (("elife-80127-supp3-v2.csv","PH1"),("elife-80127-supp6-v2.csv","PH2"),
               ("elife-80127-supp5-v2.csv","MH2"),("elife-80127-supp8-v2.csv","MH1")):
    d = pd.read_csv(os.path.join(REF,f), encoding="utf-8-sig"); d = d[d.Qvalue<0.05]
    sets[lab] = d["Gene ID"].astype(str).map(num).unique().tolist()
# 剔除细胞周期 / 管家类基因
CYCPAT = re.compile(r"(?i)^(ccn|cdk|cdc|mcm|orc|pola|prim|fen1|lig1|tyms|rrrn|dhrf|gins|pcna|aur|plk|bub|mad|top2|tpx2|kif|cenp|smc|anln|racgap|ect2|cks|tubb|tuba|espl|ndc80|nuf2|incenp|prc1|aspm|mki67|ube2|h1-|h2a|h2b|h3|h4|histone|rpl|rps|rpL|rpS|eef|act|tub)")
def filt(ids):
    out_ = []
    for x in ids:
        s = id2sym.get(x, x)
        if CYCPAT.match(str(s)): continue
        out_.append(x)
    return out_
sets_nc = {k: filt(v) for k, v in sets.items()}
for k in sets: print("%-4s marker %4d -> 剔除细胞周期/管家后 %4d" % (k, len(sets[k]), len(sets_nc[k])))
def score(A, label, groupcol, ids, key):
    A = A.copy(); A.X = A.layers["counts"].copy() if "counts" in A.layers else A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    s2i = {}
    if "sym" in A.var.columns:
        for i, s in zip(A.var_names, A.var["sym"]):
            if s is not None and pd.notna(s): s2i.setdefault(str(s), i)
    else:
        s2i = {id2sym[k2]: i for i, k2 in enumerate(A.var_names) if k2 in id2sym}
    idx = sorted({s2i[id2sym.get(x, x)] for x in ids if id2sym.get(x, x) in s2i})
    labs = [c for c in set(A.obs[groupcol].astype(str)) if str(c).isdigit()]
    rows = []
    for c in sorted(labs, key=int):
        m = (A.obs[groupcol].astype(str) == c).values
        if m.sum() < 20: continue
        rows.append({"dataset": label, "cluster": c, "n": int(m.sum()), key: float(np.asarray(A[m][:, idx].X.mean(0)).ravel().mean())})
    return pd.DataFrame(rows), len(idx)
M = sc.read_h5ad(out + r"\adata_soupx.h5ad"); M.var["sym"] = feats.set_index("gid").reindex(M.var_names)["sym"].values
E = sc.read_h5ad(out + r"\elife_adata.h5ad")
G = sc.read_h5ad(out + r"\geo_GSE306674_processed.h5ad"); G.var["sym"] = [id2sym.get(str(x), None) for x in G.var_names]
tab = {}
for lab, ids in sets_nc.items():
    r1, n1 = score(M, "我的数据", "vendor_cluster", ids, lab)
    r2, n2 = score(E, "eLife血细胞", "elife_cluster", ids, lab)
    r3, n3 = score(G, "GEO肝胰腺", "cl", ids, lab)
    tab[lab] = (r1, r2, r3, n1)
    print("\n########## %s（剔除细胞周期后，实际用 %d 个基因）" % (lab, n1))
    print("  我的数据:"); print(r1.sort_values(lab, ascending=False)[["cluster","n",lab]].round(3).to_string(index=False))
    print("  eLife（正对照，PH 应最高）:"); print(r2.sort_values(lab, ascending=False)[["cluster","n",lab]].round(3).head(5).to_string(index=False))
    print("  GEO 肝胰腺（应全低）:"); print(r3.sort_values(lab, ascending=False)[["cluster","n",lab]].round(3).head(3).to_string(index=False))
