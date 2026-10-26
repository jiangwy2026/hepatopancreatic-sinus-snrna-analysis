import os,collections
import numpy as np,pandas as pd,anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
e=ad.read_h5ad(os.path.join(OUT,"elife_adata.h5ad"))
print(e.shape, e.obs.columns.tolist())
if "counts" in e.layers: e.X=e.layers["counts"].copy()
import scanpy as sc
sc.pp.normalize_total(e,target_sum=1e4); sc.pp.log1p(e)
sym=e.var["sym"].astype(str)
def map_set(df,n=50):
    d=df.sort_values("Log2FC",ascending=False).head(n); out=[]
    for x in d["Gene Name"].astype(str):
        m=np.where(sym.str.lower().values==x.lower())[0]
        if len(m)==1 and m[0] not in out: out.append(m[0])
    return out
p1=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp3-v2.csv")),50)
p2=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp6-v2.csv")),50)
def score(idx):
    N=e.n_vars;k=len(idx);out=np.zeros(e.n_obs,dtype=np.float32)
    for st in range(0,e.n_obs,500):
        en=min(e.n_obs,st+500);X=np.asarray(e[st:en].X.todense(),dtype=np.float32)
        r=np.argsort(-X,axis=1,kind="stable");rr=np.empty_like(r);rr[np.arange(en-st)[:,None],r]=np.arange(1,N+1)
        S=rr[:,idx].sum(1);out[st:en]=1-(S-k*(k+1)/2)/(k*(N-k))
    return out
e.obs["PH1"]=score(p1);e.obs["PH2"]=score(p2)
tab=e.obs.groupby("elife_cluster").agg(n=("phase","size"),G2M=("phase",lambda x:(x=="G2M").mean()*100),PH1=("PH1","mean"),PH2=("PH2","mean")).sort_values("PH2",ascending=False)
print(tab.to_string())
