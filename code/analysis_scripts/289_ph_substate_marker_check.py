import os
import numpy as np,pandas as pd,anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
a=ad.read_h5ad(os.path.join(OUT,"step3","step3_lyo_scores.h5ad"),backed="r")
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv")); ph.columns=["cell","group","sub","u1","u2"]
sub=ph.set_index("cell")["sub"]
keep=np.isin(a.obs_names,sub.index)
A=a[keep]; s=sub.reindex(A.obs_names).values
g2v={str(g).replace("ncbi_",""):g for g in A.var_names}
genes={"PEN-3":"113808997","LOC113819543":"113819543","LOC113819837":"113819837","LOC113830563":"113830563","LOC113802679":"113802679",
"LOC113808957":"113808957","hdc":"113803948","LOC113808628":"113808628","Lbr":"113802679","LOC113800422":"113800422",
"LOC113810024":"113810024","LOC113823090":"113823090","LOC113823089":"113823089","TGM1":"113823945",
"PPAF2":"113812814","Igfbp7":"113811702","Grik2":"113823802","ccnd2":"113814652","THBS4":"113823118"}
rows=[]
for name,gid in genes.items():
    if gid not in g2v: continue
    x=np.asarray(A[:,g2v[gid]].layers["lognorm"].todense()).ravel()
    for s_ in [0,1,2,3]:
        m=s==s_
        rows.append({"gene":name,"sub":s_,"mean":x[m].mean(),"pct":(x[m]>0).mean()*100})
d=pd.DataFrame(rows)
print(d.pivot(index="gene",columns="sub",values="pct").round(1).to_string())
print()
print(d.pivot(index="gene",columns="sub",values="mean").round(2).to_string())
