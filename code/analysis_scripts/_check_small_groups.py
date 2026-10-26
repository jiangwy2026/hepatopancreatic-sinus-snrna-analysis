import io, sys, os
import pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
M = pd.read_csv(os.path.join(OUT, "markers_10groups.csv"))
for g in ("M", "MyoEpi"):
    d = M[(M.grp == g) & (M.pvals_adj < 0.05) & (M.logfoldchanges > 0.5)].copy()
    d = d.sort_values("pct_nz_group", ascending=False)
    print("== %s（按检出率排序前 10）==" % g)
    print(d.head(10)[["sym","names","pct_nz_group","pct_nz_reference","logfoldchanges"]].round(3).to_string(index=False))
