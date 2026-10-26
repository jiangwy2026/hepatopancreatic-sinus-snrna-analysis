import gzip, pandas as pd, os
f=r"D:\桌面\single cell\analysis_output\GCF_003789085.1_feature_table.txt.gz"
with gzip.open(f,'rt',encoding='utf-8',errors='replace') as fh:
 for i in range(5): print(fh.readline()[:1000])
