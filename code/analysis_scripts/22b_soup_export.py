import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; RIN = out + r"\r_input"
SOUP_MAX = 10
for s in ("blank","treated"):
    R = sp.load_npz(out + r"\_tmp_raw_%s.npz" % s).tocsr()
    tot = np.asarray(R.sum(1)).ravel()
    empty = (tot > 0) & (tot <= SOUP_MAX)
    cnt = np.asarray(R[empty].sum(0)).ravel().astype(np.float64)
    est = cnt / cnt.sum()
    pd.DataFrame({"est": est, "counts": cnt}).to_csv(RIN + r"\soup_%s.tsv" % s, sep="\t")
    print("%s: 空液滴 %d 个, 总 counts %.0f, soup 谱非零 %d 基因" % (s, empty.sum(), cnt.sum(), (cnt>0).sum()))
    del R
