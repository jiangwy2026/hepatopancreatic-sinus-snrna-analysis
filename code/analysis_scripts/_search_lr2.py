import pandas as pd
a=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
a.columns=[str(c).strip() for c in a.columns]
for term in ['Wnt','wnt','C3','complement C3','C3ar','Itga','Itgb','integrin','laminin','Lamc','Lama','Lamb','Fn1','fibronectin','Sdc2','CD36','Flt1','Fgfr','Fgf16','Pdgf','Vegf','Notch1','Delta-like','Serrate']:
 m=a[a.Symbol.astype(str).str.lower().str.contains(term.lower(),regex=False)|a.Description.astype(str).str.lower().str.contains(term.lower(),regex=False)]
 if len(m): print("\n##",term,len(m)); print(m[["GeneID","Symbol","Description"]].head(10).to_string(index=False))
