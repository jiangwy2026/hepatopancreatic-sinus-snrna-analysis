# -*- coding: utf-8 -*-
"""202_complete_table1_mapping.py - map Table 1 XP IDs using the NCBI feature table."""
import os, sys, gzip, re
import pandas as pd, numpy as np, scanpy as sc, anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
FT=os.path.join(OUT,"GCF_003789085.1_feature_table.txt.gz")
TAB=r"D:\桌面\single cell\文章所需\Table1.xlsx"
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]

def main():
    # Table1
    t=pd.read_csv(os.path.join(OUT,"Table1_rows.csv"))
    t=t.rename(columns={t.columns[0]:"section",t.columns[1]:"human_gene",t.columns[2]:"ensp",t.columns[3]:"human_cluster",t.columns[4]:"shrimp_gene",t.columns[5]:"xp_id",t.columns[6]:"haemospasia_tpm",t.columns[7]:"untreated_tpm",t.columns[8]:"identity"})
    t=t[t.xp_id.astype(str).str.startswith("XP_")].copy()
    # feature table
    ft=pd.read_csv(FT,sep="\t",dtype=str)
    ft.columns=[c.replace("# ","").strip() for c in ft.columns]
    # collect XP -> gene
    rows=[]
    for r in ft.itertuples():
        vals=[getattr(r,"product_accession",None),getattr(r,"non_redundant_refseq",None),getattr(r,"related_accession",None)]
        for v in vals:
            if isinstance(v,str) and v.startswith("XP_"):
                rows.append({"xp_id":v.split(".")[0],"gene_id":str(r.GeneID),"symbol":str(r.symbol),"name":str(r.name),"feature":str(r.feature)})
    F=pd.DataFrame(rows).drop_duplicates(["xp_id","gene_id"])
    F["xp_base"]=F.xp_id
    t["xp_base"]=t.xp_id.astype(str).str.split(".").str[0]
    M=t.merge(F[["xp_base","gene_id","symbol","name"]].drop_duplicates("xp_base"),on="xp_base",how="left")
    M.to_csv(os.path.join(OUT,"Table1_full_mapping.csv"),index=False,encoding="utf-8-sig")
    print("Table1 rows",len(M),"mapped",M.gene_id.notna().sum())
    print(M[["human_gene","human_cluster","shrimp_gene","xp_id","gene_id","symbol"]].to_string(index=False))
    # single-cell localization
    A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values; gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
    rows=[]
    for r in M[M.gene_id.notna()].itertuples():
        gid="ncbi_"+str(r.gene_id)
        if gid not in gid2i: continue
        v=A[:,gid2i[gid]].X.toarray().ravel()
        for g in ORDER:
            m=ct==g; rows.append({"human_gene":r.human_gene,"human_cluster":r.human_cluster,"shrimp_gene":r.shrimp_gene,"gene_id":gid,"symbol":r.symbol,"group":g,"mean":float(v[m].mean()),"pct":float((v[m]>0).mean()*100)})
    E=pd.DataFrame(rows); E.to_csv(os.path.join(OUT,"Table1_full_expression.csv"),index=False,encoding="utf-8-sig")
    print("expression rows",len(E),"genes",E.gene_id.nunique() if len(E) else 0)

if __name__=="__main__": main()

