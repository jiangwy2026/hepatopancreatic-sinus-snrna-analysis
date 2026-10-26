import sys, io, os, gzip, pickle
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 280)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; RIN = out + r"\r_input"
SRC = r"E:\workspace\sc_soupx\r_input"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
genes = feats["gid"].astype(str).values
sym2i = {}
for i, x in enumerate(feats["sym"].astype(str).values): sym2i.setdefault(x, i)
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
panels = {"PEN-3":"AMP","LOC113801825":"AMP(crustin)","LOC113800363":"AMP(ALF)","LOC113813598":"AMP(IB)",
          "LOC113809538":"表皮marker","A2M":"颗粒血细胞marker","LOC113823945":"增殖marker","ccnd2":"增殖marker",
          "Prss1":"F细胞marker","Try3":"F细胞marker","ACE":"R细胞marker","Oacyl":"R细胞marker","Mal-A2":"B细胞marker",
          "RASGRF1":"神经marker","ND1":"线粒体"}
all_rows = []; store = {}
for s in ("blank","treated"):
    M = scipy.io.mmread(os.path.join(SRC, "corrected_%s.mtx" % s))
    C = sp.csr_matrix(M.T).astype(np.float32)
    cells = [l.strip().lstrip("\ufeff") for l in open(os.path.join(SRC, "cells_%s.tsv" % s), encoding="utf-8") if l.strip()]
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr()
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    order = cc.sort_values("idx")["idx"].to_numpy()
    assert len(cells) == len(order) == C.shape[0], (len(cells), len(order), C.shape)
    Fsub = F[order]; clu = cc.set_index("idx").loc[order, "Cluster"]
    print("\n########## %s ##########" % s)
    print("总 counts %.0f -> %.0f (减少 %.1f%%)" % (Fsub.sum(), C.sum(), (1-C.sum()/Fsub.sum())*100))
    rows = []
    for g, tag in panels.items():
        if g not in sym2i: continue
        i = sym2i[g]
        a = np.asarray(Fsub[:, i].todense()).ravel(); b = np.asarray(C[:, i].todense()).ravel()
        for c in ["0","1","2","4","5","6","9"]:
            m = (clu == c).to_numpy()
            if m.sum() == 0: continue
            rows.append(dict(gene=g, tag=tag, cluster=c,
                             det_before=int(round((a[m]>0).mean()*100)), det_after=int(round((b[m]>0).mean()*100)),
                             mean_before=round(float(a[m].mean()),2), mean_after=round(float(b[m].mean()),2)))
    d = pd.DataFrame(rows); all_rows.append(d.assign(sample=s))
    pb = d.pivot(index=["gene","tag"], columns="cluster", values="det_before")
    pa = d.pivot(index=["gene","tag"], columns="cluster", values="det_after")
    print("--- 检出率 校正前(%) ---"); print(pb.to_string())
    print("--- 检出率 校正后(%) ---"); print(pa.to_string())
    nb = np.asarray((Fsub>0).sum(1)).ravel(); na = np.asarray((C>0).sum(1)).ravel()
    print("每细胞检出基因数中位: %.0f -> %.0f (减少 %.1f%%)" % (np.median(nb), np.median(na), (1-np.median(na)/np.median(nb))*100))
    store[s] = dict(C=C, cells=cells, clu=clu.values, order=order)
pd.concat(all_rows).to_csv(out + r"\soupx_validation.csv", index=False)
pickle.dump(store, open(out + r"\_tmp_soupx_store.pkl","wb"))
print("\nwrote soupx_validation.csv + 缓存")
