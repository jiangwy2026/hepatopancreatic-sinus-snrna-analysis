import os, re, io, sys, numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
mk = pd.read_csv(os.path.join(OUT, "blood4_markers_all.csv"))
ROWS = ["PH", "GH", "MH1", "MH2"]
HK = re.compile(r"^(Rp[LS]\d|rpl-|rps-|RpLP|eef|eIF|Act\d|Act5C|Tub|Hsp\d|hsc\d|HSC70|PRSS)", re.I)
mk["hk"] = mk["sym"].astype(str).str.match(HK).fillna(False)
mk["dPct"] = (mk.pct_nz_group - mk.pct_nz_reference) * 100
sig = mk[(~mk.is_mt) & (~mk.hk) & (mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)].copy()
sig = sig.sort_values("dPct", ascending=False)
print("每组显著 marker 数(已排除线粒体/看家):", sig.groupby("grp").size().to_dict()); print()
for g in ROWS:
    d = sig[sig.grp == g].drop_duplicates("sym").head(12)
    print("== %s 最特异的 12 个 (按 本群检出率 - 其他群检出率 排序)" % g)
    for r in d.itertuples():
        print("   %-16s 本群 %5.1f%%  其他 %5.1f%%  dPct %5.1f  log2FC %.2f  q %.1e" % (
            str(r.sym)[:16], r.pct_nz_group*100, r.pct_nz_reference*100, r.dPct, r.logfoldchanges, r.pvals_adj))
    print()
sig.to_csv(os.path.join(OUT, "blood4_markers_specificity.csv"), index=False)
print("saved blood4_markers_specificity.csv")
