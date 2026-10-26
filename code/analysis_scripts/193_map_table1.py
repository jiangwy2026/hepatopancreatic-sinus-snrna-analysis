import pandas as pd, os, re
base=r"D:\桌面\single cell\文章所需"; out=r"D:\桌面\single cell\analysis_output"
t=pd.read_excel(os.path.join(base,"Table1.xlsx"),header=0)
# clean columns manually
t=t.rename(columns={t.columns[0]:"section",t.columns[1]:"human_gene",t.columns[2]:"ensp",t.columns[3]:"human_cluster",t.columns[4]:"shrimp_gene",t.columns[5]:"xp_id",t.columns[6]:"haemospasia_tpm",t.columns[7]:"untreated_tpm",t.columns[8]:"identity"})
t=t[t["xp_id"].notna()].copy()
print("table1 rows",len(t),"unique xp",t.xp_id.nunique())
st=pd.read_excel(os.path.join(base,"Supplemental Table 5.xlsx"),sheet_name="tpm>100")
st.columns=[str(c).strip() for c in st.columns]
print("st cols",st.columns.tolist())
# protein_id / id mappings
maps=[]
for _,r in t.iterrows():
    xp=str(r.xp_id).strip()
    hit=st[st["protein_id"].astype(str).str.strip()==xp]
    if len(hit):
        rr=hit.iloc[0]
        maps.append({"human_gene":r.human_gene,"human_cluster":r.human_cluster,"shrimp_gene":r.shrimp_gene,"xp_id":xp,
                     "gene_id":rr.get("id"),"symbol":rr.get("Symbol"),"match":"XP protein_id"})
    else:
        # fallback by shrimp symbol
        hit2=st[st["Symbol"].astype(str).str.lower()==str(r.shrimp_gene).lower()]
        if len(hit2):
            rr=hit2.iloc[0]
            maps.append({"human_gene":r.human_gene,"human_cluster":r.human_cluster,"shrimp_gene":r.shrimp_gene,"xp_id":xp,
                         "gene_id":rr.get("id"),"symbol":rr.get("Symbol"),"match":"symbol"})
        else:
            maps.append({"human_gene":r.human_gene,"human_cluster":r.human_cluster,"shrimp_gene":r.shrimp_gene,"xp_id":xp,
                         "gene_id":None,"symbol":None,"match":"unmapped"})
M=pd.DataFrame(maps)
M.to_csv(os.path.join(out,"Table1_ortholog_mapping.csv"),index=False,encoding="utf-8-sig")
print(M.to_string(index=False))
print("mapped",sum(M.match!="unmapped"),"/",len(M))
