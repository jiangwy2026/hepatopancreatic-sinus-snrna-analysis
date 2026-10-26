# -*- coding: utf-8 -*-
"""Build Tables S13-S18 source data and JSON."""
from pathlib import Path
import json, math
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT=Path(r"D:\桌面\single cell"); A=ROOT/"analysis_output"; OUT=ROOT/"修订交付_2026-10-26_v2"/"source_tables"; OUT.mkdir(parents=True,exist_ok=True)
SAMPLE_MAP={"blank":"Hs_control","treated":"Hs_hemospasia"}; COND_MAP={"blank":"control","treated":"hemospasia"}; STATE_ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]; N=26778; TABLES=[]

def clean(v):
    if v is None:return None
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,(np.floating,float)):return None if not np.isfinite(v) else float(v)
    if isinstance(v,(np.bool_,)):return bool(v)
    if pd.isna(v):return None
    return v

def add(sheet,title,caption,df,footnotes=None):
    df=df.copy(); df.to_csv(OUT/f"{sheet}.csv",index=False,encoding="utf-8-sig")
    TABLES.append({"sheet":sheet,"title":title,"caption":caption,"headers":df.columns.tolist(),"rows":[[clean(x) for x in r] for r in df.itertuples(index=False,name=None)],"n_rows":len(df),"n_cols":df.shape[1],"footnotes":footnotes or []}); print(sheet,df.shape)

print("load")
a=ad.read_h5ad(A/"step3"/"step3_lyo_scores.h5ad"); traj=ad.read_h5ad(A/"step2"/"step2_integrated.h5ad"); ph=ad.read_h5ad(A/"figures_pub"/"2026-10-23_Fig2D_PH_tSNE_v2"/"PH_only_tsne.h5ad"); ln=ad.read_h5ad(A/"adata_soupx_lognorm.h5ad")
obs=a.obs.copy(); obs["sample_id"]=obs["sample"].astype(str).map(SAMPLE_MAP); obs["condition"]=obs["sample"].astype(str).map(COND_MAP)
sym_to_var={str(s):str(v) for s,v in zip(ln.var["sym"],ln.var_names)}

print("DE")
frames=[]
for state in STATE_ORDER:
    print(" ",state)
    sub=a[a.obs.cell_type.astype(str).eq(state)].copy(); sub.obs["sample"]=pd.Categorical(sub.obs["sample"].astype(str),categories=["blank","treated"])
    sc.tl.rank_genes_groups(sub,groupby="sample",groups=["treated"],reference="blank",method="wilcoxon",pts=True,use_raw=False,n_genes=sub.n_vars)
    d=sc.get.rank_genes_groups_df(sub,group="treated"); d["pct_nz_reference"]=sub.uns["rank_genes_groups"]["pts"].loc[d.names,"blank"].values; d["gene_name"]=sub.var.loc[d.names,"sym"].values; d=d[(d.pvals_adj<0.05)&(d.logfoldchanges.abs()>0.5)].copy(); d["state_label"]=state; frames.append(d)
de=pd.concat(frames,ignore_index=True); de["state_label"]=pd.Categorical(de.state_label,STATE_ORDER,ordered=True); de=de.sort_values(["state_label","pvals_adj","logfoldchanges"],ascending=[True,True,False])
s13=pd.DataFrame({"state_label":de.state_label.astype(str),"gene_id":de.names.astype(str),"gene_name":de.gene_name.fillna(de.names).astype(str),"log2FC":de.logfoldchanges.astype(float),"p_value":de.pvals.astype(float),"p_adj_BH":de.pvals_adj.astype(float),"pct_control":de.pct_nz_reference.astype(float)*100,"pct_hemospasia":de.pct_nz_group.astype(float)*100,"test_method":"Wilcoxon rank-sum; BH across genes","notes":"Nucleus-level exploratory comparison; one pooled library per condition (n=1)."})
add("S13","Table S13. Differential expression: control vs hemospasia","Within-state treated-versus-control differential expression with BH-adjusted P<0.05 and absolute log2FC>0.5. Positive log2FC denotes higher expression in hemospasia.",s13,["Each condition is one pooled library (n=1). These nucleus-level tests are exploratory and do not substitute for biological replication."])

print("trajectory bins")
phv=ph.copy(); phv.obs["sample"]=phv.obs["sample"].astype(str); phv.obs["condition"]=phv.obs["sample"].map(COND_MAP); phv.obs["DPT"]=traj.obs.dpt_hemocyte.reindex(phv.obs_names).values; phv.obs["projected_DPT"]=traj.obs.projected_pseudotime.reindex(phv.obs_names).values; phv.obs["PH_substate"]=[f"PH_s{int(x)}" for x in phv.obs.subcluster]
heat_genes=pd.read_csv(A/"Fig5_PH_dpt_heatmap_data.csv").gene.drop_duplicates().tolist(); gix=[sym_to_var[g] for g in heat_genes if g in sym_to_var]; vlook={str(v):i for i,v in enumerate(phv.var_names)}; pos=[vlook[g] for g in gix]
phv.obs["PH_gene_set_mean"]=np.asarray(phv.X[:,pos].mean(axis=1)).ravel(); ccnd2=sym_to_var.get("ccnd2"); phv.obs["ccnd2_expression"]=np.asarray(phv[:,ccnd2].X.toarray()).ravel() if ccnd2 else np.nan
rows=[]
for sample in ["blank","treated"]:
    d=phv.obs[phv.obs["sample"].eq(sample)].copy().sort_values("DPT"); d["pseudotime_bin"]=np.floor(np.arange(len(d))*20/len(d)).astype(int)+1
    for b,g in d.groupby("pseudotime_bin",sort=True):
        rows.append({"nucleus_id_or_bin":f"PH_{COND_MAP[sample]}_bin{int(b):02d}","sample_id":SAMPLE_MAP[sample],"condition":COND_MAP[sample],"state_label":"PH","PH_substate":g.PH_substate.mode().iloc[0],"DPT":g.DPT.median(),"projected_DPT":g.projected_DPT.median(),"ccnd2_expression":g.ccnd2_expression.mean(),"bin_mean_expression_gene_set":g.PH_gene_set_mean.mean(),"notes":f"n_nuclei={len(g)}; gene_set={','.join(heat_genes)}; equal-count bin per condition"})
add("S14","Table S14. Trajectory and pseudotime source data","PH pseudotime binned into 20 equal-count bins separately for control and hemospasia. ccnd2_expression is mean log-normalized expression; bin_mean_expression_gene_set is the mean expression of the 14 PH progression genes.",pd.DataFrame(rows))

print("external programs")
gs=pd.read_csv(A/"step3"/"step3_genesets.csv").set_index("set"); prog_cols={"LyO_Progenitor":"ucell_LyO_Progenitor","LyO_Lymphocyte_precursor":"ucell_LyO_Lymphocyte_precursor","LyO_Macrophage_like":"ucell_LyO_Macrophage_like","LyO_Hemocytes":"ucell_LyO_Hemocytes","LyO_ECM_producing":"ucell_LyO_ECM_producing"}
markers=pd.read_csv(OUT/"S10.csv"); top={s:set(markers[markers.state_label.eq(s)].sort_values(["p_adj_BH","log2FC"],ascending=[True,False]).head(50).gene_id.astype(str)) for s in STATE_ORDER}
subm=pd.read_csv(OUT/"S11.csv"); subtop={s:set(subm[subm.PH_substate.eq(s)].sort_values(["p_adj_BH","log2FC"],ascending=[True,False]).head(50).gene_id.astype(str)) for s in ["PH_s0","PH_s1","PH_s2","PH_s3"]}
rows=[]
for state in STATE_ORDER:
    mask=obs.cell_type.astype(str).eq(state).values
    for prog,col in prog_cols.items():
        ext=set(str(gs.loc[prog,"genes"]).split(";")); ov=sorted(ext&top[state]); exp=len(ext)*len(top[state])/N; pv=hypergeom.sf(len(ov)-1,N,len(ext),len(top[state])) if ov else 1.0; score=obs[col].astype(float)
        rows.append({"state_label":state,"program_name":prog,"mean_UCell_score":score[mask].mean(),"pct_above_global_median":(score[mask]>score.median()).mean()*100,"overlap_genes":";".join(ov),"n_overlap":len(ov),"expected_overlap":exp,"enrichment_fold":len(ov)/exp if exp else np.nan,"hypergeom_p_value":pv,"p_adj_BH":np.nan,"notes":"UCell-style rank score; overlap uses the LyO program set and the state top-50 significant markers."})
ext=set(str(gs.loc["LyO_Progenitor","genes"]).split(";")); score=phv.obs.ucell_LyO_Progenitor.astype(float)
for sub in ["PH_s0","PH_s1","PH_s2","PH_s3"]:
    mask=phv.obs.PH_substate.eq(sub).values; ov=sorted(ext&subtop[sub]); exp=len(ext)*len(subtop[sub])/N; pv=hypergeom.sf(len(ov)-1,N,len(ext),len(subtop[sub])) if ov else 1.0
    rows.append({"state_label":sub,"program_name":"LyO_Progenitor","mean_UCell_score":score[mask].mean(),"pct_above_global_median":(score[mask]>score.median()).mean()*100,"overlap_genes":";".join(ov),"n_overlap":len(ov),"expected_overlap":exp,"enrichment_fold":len(ov)/exp if exp else np.nan,"hypergeom_p_value":pv,"p_adj_BH":np.nan,"notes":"PH-substate top-50 marker overlap; descriptive external-program concordance."})
s15=pd.DataFrame(rows); s15["p_adj_BH"]=multipletests(s15.hypergeom_p_value.astype(float),method="fdr_bh")[1]
add("S15","Table S15. External progenitor program scoring","UCell-style scores and top-marker hypergeometric overlaps for external lymphoid-organ programs across Hs states, plus the LyO progenitor program across PH substates.",s15,["Scores and overlaps are descriptive and do not establish cell identity or lineage equivalence."])

print("LR communication")
lr=pd.read_csv(A/"step1"/"step1_lr_communication_sonicparanoid.csv"); sub=pd.read_csv(A/"figures_pub"/"Supplementary_niche"/"Supp11_subsampling_stability.csv"); smap={(r.source,r.target):r.top3_freq/100 for r in sub.itertuples(index=False)}
vlook={str(v).replace("ncbi_",""):v for v in a.var_names}; ctypes=a.obs.cell_type.astype(str).values; treated=(a.obs["sample"].astype(str).values=="treated"); rng=np.random.default_rng(20261026)
def expr(ids):
    ix=[vlook[str(x)] for x in ids if str(x) in vlook]; return np.asarray(a[:,ix].X.mean(axis=1)).ravel() if ix else np.zeros(a.n_obs)
rows=[]
for r in lr.itertuples(index=False):
    lig=expr(str(r.ligand_genes).split(";")); rec=expr(str(r.receptor_genes).split(";")); sm=ctypes==str(r.source); tm=ctypes==str(r.target)
    pc=lig[sm&~treated].mean()*rec[tm&~treated].mean(); ph=lig[sm&treated].mean()*rec[tm&treated].mean(); stat=math.log2((ph+1e-6)/(pc+1e-6))
    sy=treated[sm]; ty=treated[tm]; sl=lig[sm]; tr=rec[tm]; null=np.empty(500)
    for k in range(500):
        psy=rng.permutation(sy); pty=rng.permutation(ty); c=sl[~psy].mean()*tr[~pty].mean(); h=sl[psy].mean()*tr[pty].mean(); null[k]=math.log2((h+1e-6)/(c+1e-6))
    pv=(np.sum(np.abs(null)>=abs(stat))+1)/501
    rows.append({"source_state":r.source,"target_state":r.target,"ligand":r.ligand,"receptor":r.receptor,"ligand_gene_id":str(r.ligand_genes),"receptor_gene_id":str(r.receptor_genes),"communication_prob_control":pc,"communication_prob_hemospasia":ph,"log2FC":stat,"permutation_p_value":pv,"subsampling_top3_frequency":smap.get((r.source,r.target),np.nan),"orthology_support":"SonicParanoid/BLAST-supported candidate pair; see Table S17","notes":f"Group-specific probability from SoupX-corrected log-normalized expression; 500 label permutations; n=1 library per condition. Overall state-label permutation P={r.p_perm:.6g}."})
add("S16","Table S16. Ligand-receptor communication","All tested candidate source-to-target ligand-receptor combinations. Group-specific probabilities use the same log-normalized expression; P values are exploratory group-label permutation tests (500 permutations).",pd.DataFrame(rows),["subsampling_top3_frequency is a source-to-target stability metric repeated across ligand-receptor pairs sharing the same source and target.","Candidate communication does not establish functional signaling."])

print("orthology")
full=pd.read_csv(A/"step1"/"human_shrimp_niche_orthologs_full.csv").set_index("human_gene"); lookup=pd.read_csv(A/"step1"/"sonicparanoid_niche_lookup.csv").set_index("human_gene"); roles={}
for r in lr.itertuples(index=False):
    for gid in str(r.ligand_genes).split(";"): roles.setdefault((r.ligand,str(gid)),set()).add("ligand")
    for gid in str(r.receptor_genes).split(";"): roles.setdefault((r.receptor,str(gid)),set()).add("receptor")
rows=[]
for (hg,gid),rs in sorted(roles.items()):
    info=full.loc[hg] if hg in full.index else None; sp=lookup.loc[hg] if hg in lookup.index else None
    rows.append({"ligand_or_receptor":"+".join(sorted(rs)),"gene_symbol":hg,"L_vannamei_GeneID":str(gid),"ortholog_support":str(info.verdict) if info is not None else "not available","SonicParanoid_cluster":str(int(sp.group_id)) if sp is not None and pd.notna(sp.group_id) else "not detected in current SonicParanoid query set","notes":f"BLASTP best XP={info.best_xp}; product={info['product']}; identity={info.pident:.1f}%; coverage={int(info.qcovs)}%" if info is not None else "No BLASTP annotation in the candidate table"})
add("S17","Table S17. LR orthology mapping","Ligand and receptor genes used in candidate communication analysis, with BLASTP support and SonicParanoid orthogroup membership where available.",pd.DataFrame(rows),["BLASTP best hits are evidence of sequence-level orthology support, not definitive one-to-one orthology."])

print("ITGB1")
itg=[g for g in ["ncbi_113830322","ncbi_113805715"] if g in ln.var_names]; idx=[vlook[str(v).replace("ncbi_","")] for v in itg]; val=np.asarray(ln[:,idx].X.mean(axis=1)).ravel(); rows=[]
groups={s:obs.cell_type.astype(str).eq(s).values for s in STATE_ORDER}; groups.update({"Niche (TS/MH2/MyoEpi/M)":obs.cell_type.astype(str).isin(["TS","MH2","MyoEpi","M"]).values,"Hepatopancreas-like (F/R/B)":obs.cell_type.astype(str).isin(["F","R","B"]).values,"Other (PM-like)":obs.cell_type.astype(str).eq("PM-like").values,"All Hs nuclei":np.ones(a.n_obs,dtype=bool)})
for g,m in groups.items():
    x=val[m]; rows.append({"group":g,"mean_log_normalized_expression":x.mean(),"pct_positive":(x>0).mean()*100,"median_expression":np.median(x),"n_nuclei":int(m.sum()),"notes":"Mean of the two annotated ITGB1 orthologs (GeneIDs 113830322 and 113805715); positive means expression >0."})
add("S18","Table S18. ITGB1 expression summary","ITGB1 expression across Hs states and aggregate receiver groups. Expression is the per-nucleus mean of the two annotated ITGB1 orthologs.",pd.DataFrame(rows))
(OUT/"tables_S13_S18.json").write_text(json.dumps({"tables":TABLES},ensure_ascii=False,allow_nan=False),encoding="utf-8")
print("wrote",len(TABLES))
