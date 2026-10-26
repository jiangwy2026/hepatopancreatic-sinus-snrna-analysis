import sys, io, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; RIN = out + r"\r_input"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
genes = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, usecols=[0]).iloc[:,0].astype(str).values
SOUP_MAX = 10
for s in ("blank","treated"):
    R = sp.load_npz(out + r"\_tmp_raw_%s.npz" % s).tocsr()
    tot = np.asarray(R.sum(1)).ravel()
    empty = (tot > 0) & (tot <= SOUP_MAX)
    cnt = np.asarray(R[empty].sum(0)).ravel().astype(np.float64)
    est = cnt / cnt.sum()
    pd.DataFrame({"est": est, "counts": cnt}, index=genes).to_csv(RIN + r"\soup_%s.tsv" % s, sep="\t")
    print("%s: 写出 %d 行, 首行=%s, est 之和=%.4f" % (s, len(genes), genes[0], est.sum()))
    del R
print("OK")
