import pandas as pd, os
f=r"D:\桌面\single cell\analysis_output\bulk_GHRL25030101\5.ExpressionStat\all.genes.expression.annot.xls"
# inspect first lines text
for enc in ["utf-8","gb18030"]:
 try:
  with open(f,encoding=enc) as fh:
   for i in range(3): print(enc,i,fh.readline()[:1000])
  break
 except Exception as e: print("err",enc,e)
d=pd.read_csv(f,sep="\t",nrows=5)
print(d.shape); print(d.columns.tolist()); print(d.head().to_string())
