import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
sym = feats["sym"].astype(str).values
sym2idx = {}
for i, x in enumerate(sym): sym2idx.setdefault(x, i)

# 生物学上明确"某些群不该表达"的基因（高丰度、高特异）
AMP = ["PEN-3","LOC113801825","LOC113800363","LOC113813598","LOC113815775","LOC113820206","LOC113816267"]
CUT = ["LOC113809538","LOC113809539","LOC113809540","LOC113819135"]
DIG = ["Prss1","Prss2","Try3","Jon99Cii","Cpa2","LCP2"]
HP  = ["ACE","Oacyl","Mal-A2"]
hemocyte_clusters = {"0","1","2"}
cuticle_clusters  = {"3","4"}
digest_clusters   = {"5","6","7"}

cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)

for s in ("blank","treated"):
    print("\n################ %s ################" % s)
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr()
    soup = np.load(out + r"\_tmp_soup_%s.npy" % s)
    tot = np.asarray(F.sum(1)).ravel()
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    clu = cc.set_index("idx")["Cluster"]
    rows = []
    for c in sorted(set(clu), key=int):
        idx = clu[clu==c].index.to_numpy()
        if len(idx) < 30: continue
        U = tot[idx].mean(); sub = F[idx]
        obs = np.asarray(sub.mean(0)).ravel().astype(float)
        cand = []
        for gl, should_exp in ((AMP, c in hemocyte_clusters), (CUT, c in cuticle_clusters),
                               (DIG, c in digest_clusters), (HP, c in {"6","7"})):
            if should_exp: continue
            for g in gl:
                if g not in sym2idx: continue
                i = sym2idx[g]; a = soup[i]*U
                if a < 0.5: continue            # 期望太低无法估计
                cand.append(obs[i]/a)
        if len(cand) < 3:
            rows.append(dict(cluster=c, n=len(idx), UMI=round(U), n_gene=len(cand), rho=None)); continue
        rows.append(dict(cluster=c, n=len(idx), UMI=round(U), n_gene=len(cand),
                         rho=round(float(np.median(cand)),3),
                         rho_q25=round(float(np.percentile(cand,25)),3),
                         rho_q75=round(float(np.percentile(cand,75)),3)))
    d = pd.DataFrame(rows); print(d.to_string(index=False))
    v = d.rho.dropna()
    print("  可用群 %d 个 | rho 中位 %.3f | 范围 %.3f - %.3f" % (len(v), v.median(), v.min(), v.max()))
