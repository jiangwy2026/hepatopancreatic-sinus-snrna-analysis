import pandas as pd, os, re
f=r"D:\桌面\single cell\文章所需\Supplemental Table 5.xlsx"
d=pd.read_excel(f,sheet_name="tpm>100")
print(d[["id","Symbol","untreated_tpm","haemospasia_tpm","protein_id"]].head(10).to_string(index=False))
print("n",len(d),"unique id",d.id.nunique())
# clean id to ncbi_
def clean(x):
 s=str(x)
 m=re.search(r'(ncbi_\d+)',s)
 return m.group(1) if m else s
d["gene_id"]=d["id"].map(clean)
print(d.gene_id.head(10).tolist())
d[["gene_id","Symbol","untreated_tpm","haemospasia_tpm","Description","TF_family"]].to_csv(r"D:\桌面\single cell\analysis_output\prior_bulk_tpm_genes.csv",index=False,encoding="utf-8-sig")
