import pandas as pd, anndata as ad, collections, sys
sys.stdout.reconfigure(encoding="utf-8")
a=ad.read_h5ad(r"D:\桌面\single cell\analysis_output\step2\step2_integrated.h5ad", backed="r")
vn=set(map(str,a.var_names))
sym=collections.defaultdict(list)
for i,s in enumerate(a.var["sym"].astype(str)): sym[s.lower()].append(a.var_names[i])
d=pd.read_excel(r"D:\桌面\single cell\DeGene.list.xls", sheet_name="Top50.DeGene.list")
for cl in [13,4,2,7,0]:
    g=d[d.Target_Cluster==cl].sort_values("Log2FC",ascending=False).head(50)
    exact=uni=amb=none=0
    for x in g["Gene Name"].astype(str):
        if x in vn: exact+=1
        elif x.lower() in sym:
            if len(sym[x.lower()])==1: uni+=1
            else: amb+=1
        else: none+=1
    print(f"cluster {cl}: exact={exact} unique_symbol={uni} ambiguous={amb} unmapped={none}")
