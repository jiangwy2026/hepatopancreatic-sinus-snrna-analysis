import io, sys, os
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
obs = pd.read_csv(os.path.join(OUT, "all_19795_celltype.csv"), index_col=0)
sx = pd.read_hdf  # placeholder
# 厂商保留的 9,658 个细胞（用于"干净"地估计各类型阈值）
vendor_cells = set(pd.read_csv(os.path.join(r"D:\桌面\single cell",
    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
    "2.Cluster", "2.cluster", "Cells.cluster.list.xls"), sep="\t")["Cells"].astype(str))
obs["vendor_kept"] = [str(x) in vendor_cells for x in obs.index]
obs = obs.rename(columns={"pct_counts_mt": "mt", "n_genes_by_counts": "ng"})
print("%-5s %6s %8s %8s %10s %10s" % ("type", "n_all", "mt_med*", "cut*", "kept_ruleA", "kept_ruleB"))
res = []
for g in ORD:
    d = obs[obs["ct"] == g]
    clean = d[d.vendor_kept]
    med, mad = clean.mt.median(), np.median(np.abs(clean.mt - clean.mt.median())) * 1.4826
    cut = med + 3 * mad
    keepA = (d.mt < cut) & (d.ng >= 200)                     # 阈值来自干净细胞，应用到全部
    keepB = (d.mt < 20) & (d.ng >= 200)                      # 固定 20% 上限
    res.append({"ct": g, "n_all": len(d), "cut": round(cut, 2),
                "keptA": int(keepA.sum()), "keptB": int(keepB.sum())})
    print("%-5s %6d %8.2f %8.2f %10d %10d" % (g, len(d), med, cut, keepA.sum(), keepB.sum()))
R = pd.DataFrame(res)
print("\n合计: ruleA(按类型, 阈值取自干净细胞) %d | ruleB(mt<20%%) %d | 厂商 %d" % (R.keptA.sum(), R.keptB.sum(), len(vendor_cells)))

def compose(keepmask, label):
    d = obs[keepmask]
    t = pd.crosstab(d["ct"], d["sample"])
    p = t / t.sum(axis=0) * 100
    row = {"方案": label, "blank_n": int((d["sample"]=="blank").sum()), "treated_n": int((d["sample"]=="treated").sum())}
    for g in ORD:
        row[g] = "%.1f/%.1f (%.2fx)" % (p.loc[g,"blank"], p.loc[g,"treated"], p.loc[g,"treated"]/p.loc[g,"blank"])
    return row
rows = []
for g in ORD:
    d = obs[obs["ct"] == g]; clean = d[d.vendor_kept]
    med, mad = clean.mt.median(), np.median(np.abs(clean.mt - clean.mt.median())) * 1.4826
    pass
keepA = np.zeros(len(obs), bool)
for g in ORD:
    m = (obs["ct"] == g).values; d = obs[m]; clean = d[d.vendor_kept]
    cut = clean.mt.median() + 3 * np.median(np.abs(clean.mt - clean.mt.median())) * 1.4826
    keepA[m] = ((d.mt < cut) & (d.ng >= 200)).values
keepB = ((obs.mt < 20) & (obs.ng >= 200)).values
keepV = obs.vendor_kept.values
print("\n=== 三种 QC 方案下的组成（blank%/treated% (倍数)）===")
out = pd.DataFrame([compose(keepV, "厂商 mt<10%"), compose(keepA, "按类型 3MAD(干净阈值)"), compose(keepB, "全局 mt<20%")])
print(out.set_index("方案").T.to_string())
out.to_csv(os.path.join(OUT, "composition_qc_three_schemes.csv"), index=False)
