# -*- coding: utf-8 -*-
"""204_lyo_cross_organ_mapping.py
Marker-level mapping of lymphoid-organ clusters to hepatic sinus, hemolymph and hepatopancreas.
"""
import os, sys, re, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
LYO=r"D:\桌面\single cell\DeGene.list.xls"
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
NAME={0:"ECM-producing",1:"Macrophage-like 1",2:"Macrophage-like 2",3:"Unknown",4:"Lymphocyte precursor-like",5:"Capsular",6:"Unknown",7:"Hemocytes",8:"SLC-rich",9:"Unknown",10:"Unknown",11:"Epithelial",12:"Unknown",13:"Progenitor",14:"Unknown",15:"Unknown"}

def load_lyo():
 d=pd.read_excel(LYO,sheet_name="Top50.DeGene.list")
 d["Target_Cluster"]=d.Target_Cluster.astype(int)
 sig={}
 for c,g in d.groupby("Target_Cluster"):
  g=g.sort_values("Log2FC",ascending=False)
  syms=[]
  for x in g["Gene Name"].astype(str):
   if x.startswith("ncbi_"): syms.append(x)
   else: syms.append(x)
  sig[c]=list(dict.fromkeys(syms))[:50]
 return d,sig

def dataset_scores(A, group_col, sig, sym_map):
 # sym_map: synonym -> var index; group labels
 A=A.copy(); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 groups=[str(x) for x in A.obs[group_col].astype(str).values]
 uniq=[x for x in pd.unique(pd.Series(groups)) if x!="nan"]
 allg=sorted(set(g for c in sig.values() for g in c if g in sym_map))
 idx=[sym_map[g] for g in allg]
 X=np.asarray(A[:,idx].X.toarray(),dtype=np.float32)
 # z-score genes across cells
 X=(X-X.mean(0,keepdims=True))/(X.std(0,keepdims=True)+1e-6)
 # group mean per gene
 gm={g:X[np.array([x==g for x in groups])].mean(0) for g in uniq}
 rows=[]
 for c,genes in sig.items():
  gl=[g for g in genes if g in sym_map]
  if len(gl)<5: continue
  v=np.array([gm[g][[allg.index(x) for x in gl]].mean() for g in uniq])
  z=(v-v.mean())/(v.std()+1e-9)
  for g,raw,zz in zip(uniq,v,z): rows.append({"lyo_cluster":c,"target_group":g,"raw_score":float(raw),"z_score":float(zz),"n_genes":len(gl)})
 return pd.DataFrame(rows)

def main():
 d,sig=load_lyo()
 # our data
 A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.obs["group"]=A.obs.vendor_cluster.astype(str).map(GROUPS)
 symmap={}
 for i,s in enumerate(A.var["sym"].astype(str)): symmap.setdefault(s,i)
 H=dataset_scores(A,"group",sig,symmap); H["organ"]="Hepatic sinus"
 # eLife
 E=ad.read_h5ad(os.path.join(OUT,"elife_adata.h5ad"))
 emap={}
 for i,s in enumerate(E.var["sym"].astype(str)): emap.setdefault(s,i)
 EL=dataset_scores(E,"elife_cluster",sig,emap); EL["organ"]="Hemolymph"
 # GEO: var names numeric; map ncbi_ to symbols using our annotation
 G=ad.read_h5ad(os.path.join(OUT,"geo_GSE306674_processed.h5ad"))
 gid2sym={}
 for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): gid2sym.setdefault(gid,s)
 gmap={}
 for i,gid in enumerate(G.var_names.astype(str)):
  key="ncbi_"+gid if not gid.startswith("ncbi_") else gid
  s=gid2sym.get(key)
  if s: gmap.setdefault(s,i)
 GP=dataset_scores(G,"cl",sig,gmap); GP["organ"]="Hepatopancreas"
 R=pd.concat([H,EL,GP],ignore_index=True)
 R.to_csv(os.path.join(OUT,"lyo_cross_organ_mapping_scores.csv"),index=False,encoding="utf-8-sig")
 # top matches per LYO cluster and organ
 for c in sorted(sig):
  print("\nCluster",c,NAME.get(c,"?"))
  for org in ["Hepatic sinus","Hemolymph","Hepatopancreas"]:
   x=R[(R.lyo_cluster==c)&(R.organ==org)].sort_values("z_score",ascending=False).head(3)
   print(" ",org,": ",", ".join("%s(%.2f)"%(r.target_group,r.z_score) for r in x.itertuples()))

if __name__=="__main__": main()


