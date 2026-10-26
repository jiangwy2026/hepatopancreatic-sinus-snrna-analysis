import pandas as pd, os, re
base=r"D:\桌面\single cell\文章所需"; out=r"D:\桌面\single cell\analysis_output"
t=pd.read_excel(os.path.join(base,"Table1.xlsx"),header=0)
t=t.rename(columns={t.columns[0]:"section",t.columns[1]:"human_gene",t.columns[2]:"ensp",t.columns[3]:"human_cluster",t.columns[4]:"shrimp_gene",t.columns[5]:"xp_id",t.columns[6]:"haemospasia_tpm",t.columns[7]:"untreated_tpm",t.columns[8]:"identity"})
t=t[t["xp_id"].notna()].copy()
# skip repeated header row
t=t[t["xp_id"].astype(str).str.startswith("XP_")]
orth=pd.read_excel(os.path.join(base,"Supplemental Table 5.xlsx"),sheet_name="shrimp-human orthologs")
recs=[]
for s in orth["shrimp"].astype(str):
 parts=s.split("|")
 if len(parts)>=2:
  recs.append({"gene_id":parts[0],"xp_id":parts[1],"length":parts[2] if len(parts)>2 else None})
O=pd.DataFrame(recs).drop_duplicates("xp_id")
print("orth xp",len(O))
M=t.merge(O,on="xp_id",how="left")
# also map symbol from single-cell later; save now
M.to_csv(os.path.join(out,"Table1_ortholog_mapping_full.csv"),index=False,encoding="utf-8-sig")
print(M[["human_gene","human_cluster","shrimp_gene","xp_id","gene_id","length"]].to_string(index=False))
print("mapped",M.gene_id.notna().sum(),"/",len(M))
