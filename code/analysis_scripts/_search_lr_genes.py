import pandas as pd, re
a=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
a.columns=[str(c).strip() for c in a.columns]
terms=['notch','delta','serrate','dll','igf','igfbp','thbs','thrombospondin','itga','itgb','integrin','laminin','lam','fibronectin','fn1','complement','c3','c1q','cxcl','ccl','interleukin','tnf','tgf','bmp','wnt','fzd','egf','fgf','pdgf','vegf','cd36','syndecan','sdc']
for term in terms:
 m=a[a.Symbol.astype(str).str.lower().str.contains(term,regex=False)|a.Description.astype(str).str.lower().str.contains(term,regex=False)]
 if len(m):
  print("\n##",term,len(m)); print(m[["GeneID","Symbol","Description"]].head(12).to_string(index=False))
