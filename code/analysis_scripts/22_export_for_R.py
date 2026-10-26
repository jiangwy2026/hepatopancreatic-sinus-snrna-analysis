import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
RIN = out + r"\r_input"; os.makedirs(RIN, exist_ok=True)
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
gene_names = feats["gid"].astype(str).values          # R 侧用 GeneID 作行名
pd.Series(gene_names).to_csv(RIN + r"\genes.tsv", index=False, header=False)
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
meta = {}
for s in ("blank","treated"):
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr()
    soup = np.load(out + r"\_tmp_soup_%s.npy" % s)
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    idx = np.sort(cc["idx"].to_numpy())
    sub = F[idx]                      # cells (rows) x genes (cols)
    M = sp.coo_matrix(sub)
    # 写 MatrixMarket（基因 x 细胞）
    f = open(RIN + r"\toc_%s.mtx" % s, "w")
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write("%d %d %d\n" % (sub.shape[1], sub.shape[0], M.nnz))
    for r, c, v in zip(M.row, M.col, M.data):
        f.write("%d %d %d\n" % (c+1, r+1, int(v)))
    f.close()
    bars = ["%s_%s" % (s, fb[i]) for i in idx]
    pd.Series(bars).to_csv(RIN + r"\cells_%s.tsv" % s, index=False, header=False)
    clu = cc.set_index("idx").loc[idx, "Cluster"]
    pd.DataFrame({"cell": bars, "cluster": clu.values}).to_csv(RIN + r"\clusters_%s.tsv" % s, index=False, sep="\t")
    np.savetxt(RIN + r"\soup_%s.tsv" % s, np.column_stack([soup, np.asarray(F[idx].sum(0)).ravel()[0:0] if False else soup*0]),
               delimiter="\t", header="est\tcounts", comments="", fmt="%.10g")
    print("%s: 导出 %d 细胞 x %d 基因, nnz=%d" % (s, sub.shape[0], sub.shape[1], M.nnz))
print("导出完成 ->", RIN)
