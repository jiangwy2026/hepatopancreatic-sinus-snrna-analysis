# -*- coding: utf-8 -*-
"""209_second_fig_three_validations.py
1) permutation test for LyO program-organ links;
2) top20/30/50 marker sensitivity;
3) condition-specific program scores.
"""
import os, sys, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)

OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
LYO=r"D:\桌面\single cell\DeGene.list.xls"
PROGS={13:"Progenitor",2:"Macrophage-like",7:"Hemocytes",0:"ECM-producing",5:"Capsular",11:"Epithelial",8:"SLC-rich",4:"Lymphocyte precursor-like"}
ORG_DISPLAY={"Hepatic sinus":"Hepatopancreatic sinus (Hs)","Hemolymph":"Hemolymph","Hepatopancreas":"Hepatopancreas"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]

def make_group_z(A, group_col, sym_map, gene_names=None):
    A=A.copy(); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    X=A.X.tocsr()
    groups=[str(x) for x in A.obs[group_col].astype(str).values]
    gs=[]
    for g in pd.unique(pd.Series(groups)):
        if g in ("nan","None"): continue
        m=np.array(groups)==g
        if m.sum()>=20: gs.append(g)
    M=np.vstack([np.asarray(X[np.array(groups)==g].mean(0)).ravel() for g in gs])
    mu=np.asarray(X.mean(0)).ravel()
    ex2=np.asarray(X.multiply(X).mean(0)).ravel()
    sd=np.sqrt(np.maximum(ex2-mu*mu,1e-8))
    Z=(M-mu)/sd
    names = gene_names if gene_names is not None else (A.var["sym"].astype(str).tolist() if "sym" in A.var.columns else A.var_names.astype(str).tolist())
    return gs, Z, names, mu, sd

def load_all():
    H=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); H.obs["grp_final"]=H.obs.vendor_cluster.astype(str).map({"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"})
    hs_sym={}
    for i,s in enumerate(H.var["sym"].astype(str)): hs_sym.setdefault(s,i)
    hs_gid2sym={str(g):str(s) for g,s in zip(H.var_names.astype(str),H.var["sym"].astype(str))}
    hg,hz,hnames,hmu,hsd=make_group_z(H,"grp_final",hs_sym)
    E=ad.read_h5ad(os.path.join(OUT,"elife_adata.h5ad")); esym={}
    for i,s in enumerate(E.var["sym"].astype(str)): esym.setdefault(s,i)
    eg,ez,enames,emu,esd=make_group_z(E,"elife_cluster",esym)
    G=ad.read_h5ad(os.path.join(OUT,"geo_GSE306674_processed.h5ad")); gsym={}
    for i,gid in enumerate(G.var_names.astype(str)):
        key="ncbi_"+gid if not str(gid).startswith("ncbi_") else str(gid)
        s=hs_gid2sym.get(key)
        if s: gsym.setdefault(s,i)
    hp_names=[hs_gid2sym.get(("ncbi_"+str(x)) if not str(x).startswith("ncbi_") else str(x), str(x)) for x in G.var_names]
    gg,gz,gnames,gmu,gsd=make_group_z(G,"cl",gsym,gene_names=hp_names)
    return H,E,G,(hg,hz,hnames,hmu,hsd),(eg,ez,enames,emu,esd),(gg,gz,gnames,gmu,gsd),hs_gid2sym

def map_set(genes, name2idx, hs_gid2sym):
    out=[]
    for g in genes:
        s=str(g)
        if s.startswith("ncbi_"): s=hs_gid2sym.get(s,s)
        if s in name2idx: out.append(name2idx[s])
    return sorted(set(out))

def score_set(Z, idx):
    v=Z[:,idx].mean(1)
    z=(v-v.mean())/(v.std()+1e-9)
    return z

def main():
    H,E,G,HD,ED,GD,gid2sym=load_all()
    lyo=pd.read_excel(LYO,sheet_name="Top50.DeGene.list"); lyo["Target_Cluster"]=lyo.Target_Cluster.astype(int)
    marker={c:list(dict.fromkeys(lyo[lyo.Target_Cluster==c].sort_values("Log2FC",ascending=False)["Gene Name"].astype(str).tolist())) for c in PROGS}
    # datasets: name, groups, Z, names, mu
    DS=[("Hepatic sinus",*HD),("Hemolymph",*ED),("Hepatopancreas",*GD)]
    # ---- 1. permutation test (top50) ----
    rng=np.random.default_rng(42); perm_rows=[]
    for org,gs,Z,names,mu,sd in DS:
        name2idx={s:i for i,s in enumerate(names)}
        bins=pd.qcut(mu,10,labels=False,duplicates="drop") if len(mu)>100 else pd.qcut(mu,5,labels=False,duplicates="drop")
        bin_groups={b:np.where(bins==b)[0] for b in np.unique(bins)}
        for c,p in PROGS.items():
            idx=map_set(marker[c],name2idx,gid2sym)
            if len(idx)<5:
                perm_rows.append({"program":p,"organ":org,"n_genes":len(idx),"best_group":None,"score":np.nan,"p_perm":np.nan}); continue
            obs=score_set(Z,idx); bi=int(np.argmax(obs)); obs_score=float(obs[bi])
            null=[]
            for _ in range(1000):
                chosen=[]
                for j in idx:
                    b=bins[j]; pool=bin_groups.get(b,[])
                    if len(pool)==0: pool=np.arange(len(names))
                    chosen.append(int(rng.choice(pool)))
                null.append(float(np.max(score_set(Z,sorted(set(chosen))))))
            pval=(1+sum(x>=obs_score for x in null))/(1+len(null))
            perm_rows.append({"program":p,"organ":org,"n_genes":len(idx),"best_group":gs[bi],"score":obs_score,"p_perm":pval})
    P=pd.DataFrame(perm_rows); P["q_perm"]=multipletests(P.p_perm.fillna(1),method="fdr_bh")[1]
    P.to_csv(os.path.join(OUT,"Fig6_permutation_results.csv"),index=False,encoding="utf-8-sig")
    print("\nPermutation:\n",P.round(4).to_string(index=False))
    # ---- 2. marker sensitivity ----
    sens=[]
    for org,gs,Z,names,mu,sd in DS:
        name2idx={s:i for i,s in enumerate(names)}
        for n in (20,30,50):
            for c,p in PROGS.items():
                idx=map_set(marker[c][:n],name2idx,gid2sym)
                if len(idx)<5: sens.append({"program":p,"organ":org,"topN":n,"n_genes":len(idx),"best_group":None,"score":np.nan}); continue
                z=score_set(Z,idx); bi=int(np.argmax(z)); sens.append({"program":p,"organ":org,"topN":n,"n_genes":len(idx),"best_group":gs[bi],"score":float(z[bi])})
    S=pd.DataFrame(sens)
    S.to_csv(os.path.join(OUT,"Fig6_marker_sensitivity.csv"),index=False,encoding="utf-8-sig")
    print("\nSensitivity summary:\n",S.round(3).to_string(index=False))
    # ---- 3. condition-specific scores ----
    cond_rows=[]
    # Hs marked already; need score per cell using top50 markers
    def condition_scores(A, group_col, cond_map, label, gene_names):
        A=A.copy(); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
        X=A.X.tocsc(); names=gene_names; name2idx={s:i for i,s in enumerate(names)}
        zcell={}
        for c,p in PROGS.items():
            idx=map_set(marker[c],name2idx,gid2sym)
            if len(idx)<5: continue
            sub=X[:,idx].toarray(); mu=sub.mean(0); sd=sub.std(0)+1e-6
            zcell[p]=((sub-mu)/sd).mean(1)
            cond=np.array([cond_map(str(x)) for x in A.obs[group_col].astype(str).values])
            vals=[]
            for g in sorted(set(cond)):
                m=cond==g; vals.append(zcell[p][m])
                cond_rows.append({"dataset":label,"program":p,"condition":g,"n":int(m.sum()),"mean_score":float(zcell[p][m].mean())})
            if len(vals)>=2:
                try: kw=kruskal(*vals).pvalue
                except Exception: kw=np.nan
                for g in sorted(set(cond)):
                    cond_rows.append({"dataset":label,"program":p,"condition":g+"__p_kruskal","n":len(A.obs),"mean_score":float(kw)})
    # Hs
    H.obs["cond"]=H.obs["group"].astype(str)
    condition_scores(H,"cond",lambda x:x,"Hepatopancreatic sinus (Hs)", H.var["sym"].astype(str).tolist())
    # Hemolymph
    E.obs["cond"]=E.obs["sample"].astype(str)
    condition_scores(E,"cond",lambda x:x,"Hemolymph", E.var["sym"].astype(str).tolist())
    # HP: H/C/V
    hp_cond={"H1":"Healthy","H2":"Healthy","H3":"Healthy","C1":"Control","C2":"Control","V1":"Vibrio","V2":"Vibrio","V3":"Vibrio"}
    G.obs["cond"]=G.obs["sample"].astype(str).map(hp_cond)
    condition_scores(G,"cond",lambda x:hp_cond.get(x,x),"Hepatopancreas", [gid2sym.get("ncbi_"+str(x), gid2sym.get(str(x), str(x))) for x in G.var_names])
    C=pd.DataFrame(cond_rows); C.to_csv(os.path.join(OUT,"Fig6_condition_program_scores.csv"),index=False,encoding="utf-8-sig")
    print("\nCondition rows:",len(C))

if __name__=="__main__": main()
