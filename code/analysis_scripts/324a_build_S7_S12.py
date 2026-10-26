# -*- coding: utf-8 -*-
"""Build Tables S7-S12 source data and JSON."""
from pathlib import Path
import json, math
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import confint_proportions_2indep

ROOT=Path(r"D:\桌面\single cell"); A=ROOT/"analysis_output"; OUT=ROOT/"修订交付_2026-10-26_v2"/"source_tables"; OUT.mkdir(parents=True,exist_ok=True)
SAMPLE_MAP={"blank":"Hs_control","treated":"Hs_hemospasia"}; COND_MAP={"blank":"control","treated":"hemospasia"}
STATE_ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
TABLES=[]

def clean(v):
    if v is None: return None
    if isinstance(v,np.integer): return int(v)
    if isinstance(v,(np.floating,float)): return None if not np.isfinite(v) else float(v)
    if isinstance(v,(np.bool_,)): return bool(v)
    if pd.isna(v): return None
    return v

def add(sheet,title,caption,df,footnotes=None):
    df=df.copy(); df.to_csv(OUT/f"{sheet}.csv",index=False,encoding="utf-8-sig")
    TABLES.append({"sheet":sheet,"title":title,"caption":caption,"headers":df.columns.tolist(),"rows":[[clean(x) for x in r] for r in df.itertuples(index=False,name=None)],"n_rows":len(df),"n_cols":df.shape[1],"footnotes":footnotes or []})
    print(sheet,df.shape)

print("load")
adata=ad.read_h5ad(A/"step3"/"step3_lyo_scores.h5ad"); traj=ad.read_h5ad(A/"step2"/"step2_integrated.h5ad"); ph=ad.read_h5ad(A/"figures_pub"/"2026-10-23_Fig2D_PH_tSNE_v2"/"PH_only_tsne.h5ad"); lognorm=ad.read_h5ad(A/"adata_soupx_lognorm.h5ad")
assert adata.obs_names.equals(traj.obs_names) and adata.obs_names.equals(lognorm.obs_names)
obs=adata.obs.copy(); obs["sample_id"]=obs["sample"].astype(str).map(SAMPLE_MAP); obs["condition"]=obs["sample"].astype(str).map(COND_MAP); obs["cellranger_barcode"]=[x.split("_",1)[1] if "_" in x else x for x in obs.index.astype(str)]
s_genes=["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"]
g2m_genes=["CCNB3","AURKA","TOP2A","TPX2","KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","ASPM","Pka-C1","Patronin"]
sym_to_var={str(s):str(v) for s,v in zip(lognorm.var["sym"],lognorm.var_names)}
sc.tl.score_genes_cell_cycle(lognorm,s_genes=[sym_to_var[x] for x in s_genes if x in sym_to_var],g2m_genes=[sym_to_var[x] for x in g2m_genes if x in sym_to_var],random_state=0)
phase=lognorm.obs["phase"].astype(str).copy()
parts=[]
for sample in ["blank","treated"]:
    p=ROOT/"ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"/"ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"/"2.Cluster"/"1.QC"/"doubletfinder"/sample/f"DF.classify.{sample}.xls"
    d=pd.read_csv(p,sep="\t").rename(columns={"Cells":"nucleus_id","pANN":"doublet_score","classifications":"doublet_call"}); d["sample"]=sample; parts.append(d)
doublets=pd.concat(parts,ignore_index=True).set_index("nucleus_id").reindex(obs.index); assert doublets["doublet_score"].notna().all()
ct=pd.read_csv(A/"cross_organ_progenitor"/"cross_organ_percell_scores.csv.gz"); cross=ct[ct.target.eq("Hs PH")].set_index("cell").reindex(obs.index); assert cross["cycle_ucell"].notna().all()
ph_sub=pd.Series("",index=obs.index,dtype=object); ph_sub.loc[ph.obs.index]=[f"PH_s{int(x)}" for x in ph.obs["subcluster"]]

s7=pd.DataFrame([
{"sample_id":"Hs_control","library_id":"ADRL26080060_blank","condition":"control","timepoint":"0 h","species":"Litopenaeus vannamei","taxid":"6689","tissue":"Hepatopancreatic sinus (Hs)","n_shrimp_pooled":10,"nuclei_isolation_protocol":"NP-40 lysis; discontinuous 15%/22%/40% iodixanol gradient; 22%-40% interface","library_kit":"10x Genomics Chromium Next GEM Single Cell 3' HT Kit v3.1 (PN-1000348); Chip G (PN-1000120); Dual Index TT Set A (PN-3000431)","sequencing_platform":"Illumina NovaSeq X Plus","read_length":"PE150","reference_genome":"GCF_003789085.1 (ASM378908v1)","accession_id":"PRJNA1545411","notes":"One pooled library for this condition (library-level n=1); source label blank."},
{"sample_id":"Hs_hemospasia","library_id":"ADRL26080060_treated","condition":"hemospasia","timepoint":"6 h post-hemospasia","species":"Litopenaeus vannamei","taxid":"6689","tissue":"Hepatopancreatic sinus (Hs)","n_shrimp_pooled":10,"nuclei_isolation_protocol":"NP-40 lysis; discontinuous 15%/22%/40% iodixanol gradient; 22%-40% interface","library_kit":"10x Genomics Chromium Next GEM Single Cell 3' HT Kit v3.1 (PN-1000348); Chip G (PN-1000120); Dual Index TT Set A (PN-3000431)","sequencing_platform":"Illumina NovaSeq X Plus","read_length":"PE150","reference_genome":"GCF_003789085.1 (ASM378908v1)","accession_id":"PRJNA1545411","notes":"One pooled library for this condition (library-level n=1); source label treated."}])
add("S7","Table S7. Sample metadata","One row per snRNA-seq library. Each condition was represented by one pooled library prepared from 10 shrimp.",s7,["Library-level replication is n=1 per condition; no between-library variance can be estimated."])

s8=pd.DataFrame({"nucleus_id":obs.index.astype(str),"sample_id":obs.sample_id.values,"condition":obs.condition.values,"cellranger_barcode":obs.cellranger_barcode.values,"n_genes":obs.n_genes_by_counts.astype(int).values,"n_UMIs":obs.total_counts.astype(int).values,"mito_UMI_frac":(obs.pct_counts_mt.astype(float)/100).values,"doublet_score":doublets.doublet_score.astype(float).values,"leiden_cluster":obs.vendor_cluster.astype(int).values,"state_label":obs.cell_type.astype(str).values,"PH_substate":ph_sub.values,"cell_cycle_phase":phase.reindex(obs.index).astype(str).values,"DPT":traj.obs.dpt_hemocyte.reindex(obs.index).values,"projected_DPT":traj.obs.projected_pseudotime.reindex(obs.index).values,"UCell_PH1":cross.cycle_ucell.values,"UCell_PH2":cross.progenitor_ucell.values,"LyO_progenitor_UCell":obs.ucell_LyO_Progenitor.values,"identity_residual_z":cross.adjusted_progenitor_z.values,"notes":np.where(doublets.doublet_call.eq("Singlet"),"Retained Singlet; doublet_score is DoubletFinder pANN","REVIEW: retained despite DoubletFinder call")})
add("S8","Table S8. Nucleus-level metadata","One row per retained nucleus (n=9,658). DPT is defined for the PH-GH trajectory only; projected_DPT is available for other states. leiden_cluster contains the delivered vendor cluster ID because vendor-derived state labels were retained.",s8,["UCell_PH1 and UCell_PH2 are eLife-derived cell-cycle and progenitor-identity program scores.","All retained nuclei were classified as Singlets by the delivered DoubletFinder output."])

base=ROOT/"ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"/"ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
fs=pd.read_csv(base/"2.Cluster"/"1.QC"/"Filter.stat.xls",sep="\t").set_index("Samples"); seq=pd.read_csv(base/"1.Expression"/"samples.sequence.stat.xls",sep="\t").set_index("Sample"); align=pd.read_csv(base/"1.Expression"/"samples.align.stat.xls",sep="\t").set_index("Sample"); rho={"blank":0.314,"treated":0.356}
rows=[]
for sample in ["blank","treated"]:
    d=doublets[doublets["sample"].eq(sample)]
    rows.append({"sample_id":SAMPLE_MAP[sample],"raw_barcodes":int(fs.loc[sample,"before_filter_num"]),"retained_nuclei":int(fs.loc[sample,"after_filter_num"]),"median_genes_per_nucleus":float(fs.loc[sample,"after_filter_median_genes_per_cell"]),"median_UMIs_per_nucleus":float(fs.loc[sample,"after_filter_median_UMI_per_cell"]),"mito_UMI_frac_median":float(fs.loc[sample,"after_filter_median_MT_per_cell"])/100,"estimated_doublet_rate":float((d.doublet_call=="Doublet").mean()),"SoupX_contamination_frac":rho[sample],"reads_mapped_pct":float(str(align.loc[sample,"Reads Mapped Confidently to Transcriptome"]).replace("%","")),"Q30_pct":float(str(seq.loc[sample,"Q30 Bases in RNA Read"]).replace("%","")),"sequencing_saturation":float(str(seq.loc[sample,"Sequencing Saturation"]).replace("%","")),"notes":f"Doublet rate = DoubletFinder calls / Cell Ranger barcodes. Q30 is RNA-read Q30. Genome mapping was {align.loc[sample,'Reads Mapped Confidently to Genome']}; transcriptome mapping is reported in reads_mapped_pct."})
add("S9","Table S9. QC and sequencing summary","One row per library. Percent values use a 0-100 scale except mito_UMI_frac_median and SoupX_contamination_frac, which are fractions.",pd.DataFrame(rows))

m=pd.read_csv(A/"markers_10groups.csv"); m=m[(m.pvals_adj<0.05)&(m.logfoldchanges>0.5)].copy(); m["grp"]=pd.Categorical(m.grp,STATE_ORDER,ordered=True); m=m.sort_values(["grp","pvals_adj","logfoldchanges"],ascending=[True,True,False])
s10=pd.DataFrame({"state_label":m.grp.astype(str),"gene_id":m.names.astype(str),"gene_name":m.sym.fillna(m.names).astype(str),"log2FC":m.logfoldchanges.astype(float),"p_value":m.pvals.astype(float),"p_adj_BH":m.pvals_adj.astype(float),"pct_in_state":m.pct_nz_group.astype(float)*100,"pct_out_state":m.pct_nz_reference.astype(float)*100,"test_method":"Wilcoxon rank-sum; BH across genes"})
add("S10","Table S10. Marker genes for Hs states","Positive one-vs-rest markers with BH-adjusted P<0.05 and log2FC>0.5.",s10)

print("PH markers")
po=ph.copy(); po.obs["subcluster"]=pd.Categorical(po.obs["subcluster"].astype(int),categories=[0,1,2,3]); sc.tl.rank_genes_groups(po,groupby="subcluster",method="wilcoxon",pts=True,use_raw=False,n_genes=po.n_vars); pm=sc.get.rank_genes_groups_df(po,group=None); pm=pm[(pm.pvals_adj<0.05)&(pm.logfoldchanges>0.5)].copy(); pm["PH_substate"]=pm.group.astype(int).map(lambda x:f"PH_s{int(x)}"); pm["gene_name"]=po.var.loc[pm.names,"sym"].values; pm["PH_substate"]=pd.Categorical(pm.PH_substate,["PH_s0","PH_s1","PH_s2","PH_s3"],ordered=True); pm=pm.sort_values(["PH_substate","pvals_adj","logfoldchanges"],ascending=[True,True,False])
s11=pd.DataFrame({"PH_substate":pm.PH_substate.astype(str),"gene_id":pm.names.astype(str),"gene_name":pm.gene_name.fillna(pm.names).astype(str),"log2FC":pm.logfoldchanges.astype(float),"p_value":pm.pvals.astype(float),"p_adj_BH":pm.pvals_adj.astype(float),"pct_in_substate":pm.pct_nz_group.astype(float)*100,"pct_out_substate":pm.pct_nz_reference.astype(float)*100})
add("S11","Table S11. PH substate marker genes","Positive one-vs-rest PH-substate markers with BH-adjusted P<0.05 and log2FC>0.5.",s11)

comp=pd.read_csv(A/"composition_new_names_10groups.csv").set_index("group").reindex(STATE_ORDER).reset_index(); n0=int(comp.n_blank.sum()); n1=int(comp.n_treated.sum()); rows=[]
for r in comp.itertuples(index=False):
    a,b=int(r.n_treated),int(r.n_blank); p1,p2=a/n1,b/n0; pp=(a+b)/(n0+n1); se=math.sqrt(pp*(1-pp)*(1/n1+1/n0)); z=(p1-p2)/se if se else np.nan; pv=2*norm.sf(abs(z)) if np.isfinite(z) else np.nan; lo,hi=confint_proportions_2indep(a,n1,b,n0,method="newcomb",compare="diff")
    rows.append({"state_label":r.group,"control_n":b,"control_pct":p2*100,"hemospasia_n":a,"hemospasia_pct":p1*100,"fold_change":p1/p2 if p2 else np.nan,"log2FC":math.log2(p1/p2) if p1 and p2 else np.nan,"z_statistic":z,"p_value":pv,"p_adj_BH":np.nan,"CI_low":lo,"CI_high":hi,"notes":"Exploratory nucleus-level comparison; one pooled library per condition (n=1)."})
s12=pd.DataFrame(rows); s12["p_adj_BH"]=multipletests(s12.p_value.astype(float),method="fdr_bh")[1]
add("S12","Table S12. Cell-state composition statistics source data","Counts, percentages, fold change, two-proportion z statistic and 95% Newcombe confidence interval for the difference in proportions (hemospasia minus control). CI limits are on the proportion scale.",s12,["EXPLORATORY ONLY: each condition contains n=1 pooled library. Nucleus-level P values and FDR values do not estimate biological replicate variance and must not be interpreted as confirmatory group-level replication."])
(OUT/"tables_S7_S12.json").write_text(json.dumps({"tables":TABLES},ensure_ascii=False,allow_nan=False),encoding="utf-8")
print("wrote",len(TABLES))
