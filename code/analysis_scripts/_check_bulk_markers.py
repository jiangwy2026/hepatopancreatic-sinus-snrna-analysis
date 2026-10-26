import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
M=pd.read_csv(os.path.join(O,"top_markers_10groups.csv"))
B=pd.read_csv(os.path.join(O,"prior_bulk_tpm_genes.csv"))
print("M",M.shape,"B",B.shape)
for g in ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]:
 ids=M[(M.grp==g)&(M["rank"]<=10)].gid.astype(str).tolist()
 hit=[x for x in ids if x in set(B.gene_id.astype(str))]
 print(g,"ids",len(ids),"hit",len(hit),hit[:5])
print(B[B.gene_id.astype(str).str.contains("113823945|113814652|113824160|113808997")][["gene_id","Symbol","untreated_tpm","haemospasia_tpm"]])
