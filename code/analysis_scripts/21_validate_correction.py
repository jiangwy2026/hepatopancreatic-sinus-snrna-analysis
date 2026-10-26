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
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
panel = {"PEN-3":"AMP","LOC113801825":"AMP(crustin)","LOC113800363":"AMP(ALF)","LOC113813598":"AMP(IB)",
         "A2M":"颗粒血细胞marker","LOC113823945":"增殖marker","ccnd2":"增殖marker",
         "Prss1":"F细胞marker","Try3":"F细胞marker","ACE":"R细胞marker","Oacyl":"R细胞marker","Mal-A2":"B细胞marker",
         "LOC113809538":"表皮marker","RASGRF1":"神经marker","LOC113803948":"神经内分泌marker"}
res = {}
for s in ("blank","treated"):
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr(); C = sp.load_npz(out + r"\_tmp_corrected_%s.npz" % s).tocsr()
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    idx = cc["idx"].to_numpy(); clu = cc.set_index("idx")["Cluster"]
    print("\n################ %s ################" % s)
    print("=== 校正前后：检出率(%% ) / 平均counts ===")
    print("%-22s %-8s %s" % ("基因", "群", "  校正前 -> 校正后"))
    rows = []
    for g, tag in panel.items():
        if g not in sym2idx: continue
        i = sym2idx[g]
        a = np.asarray(F[idx][:, i].todense()).ravel(); b = np.asarray(C[idx][:, i].todense()).ravel()
        for c in sorted(set(clu), key=int):
            m = (clu == c).to_numpy()
            r = dict(gene=g, tag=tag, cluster=c, det_before=(a[m]>0).mean()*100, det_after=(b[m]>0).mean()*100,
                     mean_before=a[m].mean(), mean_after=b[m].mean())
            rows.append(r)
    d = pd.DataFrame(rows)
    piv = d.pivot_table(index=["gene","tag"], columns="cluster", values=["det_before","det_after"], aggfunc="first")
    show = d[d.cluster.isin(["0","1","4","5","6","9"])]
    print(show.pivot_table(index=["gene","tag"], columns="cluster", values=["det_before","det_after"]).round(0).to_string())
    d.to_csv(out + r"\correction_validation_%s.csv" % s, index=False)
    nb = np.asarray((F[idx]>0).sum(1)).ravel(); na = np.asarray((C[idx]>0).sum(1)).ravel()
    print("\n每细胞检出基因数: 中位 %.0f -> %.0f (减少 %.1f%%)" % (np.median(nb), np.median(na), (1-np.median(na)/np.median(nb))*100))
    print("每细胞总 counts: 中位 %.0f -> %.0f" % (np.median(np.asarray(F[idx].sum(1)).ravel()), np.median(np.asarray(C[idx].sum(1)).ravel())))
