import pandas as pd, re
p=r"D:\桌面\single cell\all.genes.marker.stat.xls.xls"
a=pd.read_csv(p,sep="\t",encoding="gb18030",low_memory=False)
a.columns=[str(c).strip() for c in a.columns]
terms=["eIF4","EIF4","ARHGAP12","dynein","eEF1","EEF1","RPS3A","RPS9","filamin","FLNA","parvin","PARVB","SLC6A12","BGT1","CHDH","BHMT","NFAT5","SOD2","SOD","ferritin","FTH1","GLUL","glutamine synthetase","betaine","choline","homocysteine"]
for term in terms:
 m=a[a.apply(lambda r: term.lower() in (str(r.get("Symbol",""))+" "+str(r.get("Description",""))).lower(),axis=1)]
 if len(m):
  print("\n##",term,len(m)); print(m[["GeneID","Symbol","Description","TF_family"]].head(15).to_string(index=False))
