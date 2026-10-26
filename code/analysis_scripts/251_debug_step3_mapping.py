import pandas as pd, anndata as ad, numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")
a=ad.read_h5ad(r"D:\桌面\single cell\analysis_output\step2\step2_integrated.h5ad", backed="r")
d=pd.read_excel(r"D:\桌面\single cell\DeGene.list.xls", sheet_name="Top50.DeGene.list")
g=d[d.Target_Cluster==13].sort_values("Log2FC",ascending=False).head(50)
sym={}
for i,s in enumerate(a.var["sym"].astype(str)): sym.setdefault(s.lower(), a.var_names[i])
vn=set(map(str,a.var_names))
mapped=[]
for x in g["Gene Name"].astype(str):
    v=x if x in vn else sym.get(x.lower())
    mapped.append((x,v))
print("mapped:"); [print("  ",x,"->",v) for x,v in mapped]
mk=pd.read_csv(r"D:\桌面\single cell\analysis_output\markers_10groups.csv")
ph=set(mk[mk.grp=="PH"].sort_values("logfoldchanges",ascending=False).head(50)["names"].astype(str))
A=set(v for x,v in mapped if v)
print("\nPH marker sample:",sorted(ph)[:10])
print("A size",len(A),"PH size",len(ph),"intersection",sorted(A&ph))
print("PPAF2 mapping in A:",[ (x,v) for x,v in mapped if x.lower()=="ppaf2"])
print("PPAF2 in PH:",[x for x in ph if "113812814" in x])
