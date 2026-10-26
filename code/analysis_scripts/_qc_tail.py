import gzip, os, numpy as np, scipy.io as sio
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
EXPR = os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions")
for s in ("blank", "treated"):
    m = sio.mmread(os.path.join(EXPR, s, "matrix.mtx.gz"))
    tot = np.asarray(m.sum(axis=0)).ravel()
    ng = np.asarray((m > 0).sum(axis=0)).ravel()
    print("%s 厂商 filtered 矩阵: %d 细胞 | UMI 中位 %d, 均值 %.0f, p99 %d, 最大 %d" % (
        s, m.shape[1], np.median(tot), tot.mean(), np.percentile(tot, 99), tot.max()))
    # 高 UMI 尾部
    for q in (99.9, 99.99):
        print("     p%.2f = %d" % (q, np.percentile(tot, q)))
    print("     n_genes 中位 %d, 最大 %d" % (np.median(ng), ng.max()))
