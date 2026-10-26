import pandas as pd, os, numpy as np
O=r"D:\桌面\single cell\analysis_output"
S=pd.read_csv(os.path.join(O,"lyo_cross_organ_mapping_scores.csv"))
PROGRAMS={13:"Progenitor",2:"Macrophage-like",7:"Hemocytes",0:"ECM-producing",5:"Capsular",11:"Epithelial",8:"SLC-rich",4:"Lymphocyte precursor-like"}
rows=[]
for c,n in PROGRAMS.items():
 for org in ["Hepatic sinus","Hemolymph","Hepatopancreas"]:
  x=S[(S.lyo_cluster==c)&(S.organ==org)].sort_values("z_score",ascending=False)
  if len(x):
   r=x.iloc[0]; rows.append({"program":n,"organ":org,"best_cluster":r.target_group,"raw_score":r.raw_score,"z_score":r.z_score,"n_genes":r.n_genes})
  else:
   rows.append({"program":n,"organ":org,"best_cluster":None,"raw_score":np.nan,"z_score":np.nan,"n_genes":0})
R=pd.DataFrame(rows)
print(R.round(3).to_string(index=False))
R.to_csv(os.path.join(O,"Fig6_program_organ_summary.csv"),index=False,encoding="utf-8-sig")
