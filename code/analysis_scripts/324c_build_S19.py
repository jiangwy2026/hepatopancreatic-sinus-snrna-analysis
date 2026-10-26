# -*- coding: utf-8 -*-
"""Build Table S19 figure source data (one sheet per figure)."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd

ROOT=Path(r"D:\桌面\single cell"); A=ROOT/"analysis_output"; OUT=ROOT/"修订交付_2026-10-26_v2"/"source_tables"; OUT.mkdir(parents=True,exist_ok=True)
SAMPLE_MAP={"blank":"Hs_control","treated":"Hs_hemospasia"}; TABLES=[]
COLS=["figure_panel","x_variable","y_variable","group","value","statistic","p_value","p_adj","n","notes"]

def clean(v):
    if v is None:return None
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,(np.floating,float)):return None if not np.isfinite(v) else float(v)
    if isinstance(v,(np.bool_,)):return bool(v)
    if pd.isna(v):return None
    return v

def add(sheet,title,caption,df):
    df=df[COLS].copy(); df.to_csv(OUT/f"{sheet}.csv",index=False,encoding="utf-8-sig")
    TABLES.append({"sheet":sheet,"title":title,"caption":caption,"headers":COLS,"rows":[[clean(x) for x in r] for r in df.itertuples(index=False,name=None)],"n_rows":len(df),"n_cols":len(COLS),"footnotes":[]}); print(sheet,df.shape)

# Figure 1: BrdU/FACS source data.
st1=pd.read_excel(ROOT/"supplement table"/"Supplemental Table 1.xlsx",sheet_name="Summary ")
rows=[]
for r in st1.itertuples(index=False,name=None):
    tissue=str(r[0]); c=[float(x) for x in r[1:5]]; h=[float(x) for x in r[5:9]]; p=float(r[10])
    for grp,vals in [("control",c),("hemospasia",h)]:
        rows.append({"figure_panel":"Fig1E","x_variable":"tissue","y_variable":"BrdU-positive cells (%)","group":f"{tissue}|{grp}","value":float(np.mean(vals)),"statistic":"mean of four biological replicates","p_value":p,"p_adj":None,"n":len(vals),"notes":"Source: Supplemental Table 1, Summary sheet; n=4 independent batches."})
add("S19_Fig1","Table S19-1. Figure 1 source data","BrdU/FACS source summary used for Figure 1E.",pd.DataFrame(rows))

# Figure 2: compartment comparison summary.
rows=[]
for label,value,n in [("Hsp-up metabolites",244,6),("Hsp-down metabolites",39,6),("Cp-only proteins",634,6),("Hsp-only proteins",163,6),("Shared proteins",67,6)]:
    rows.append({"figure_panel":"Fig2A","x_variable":"plasma compartment comparison","y_variable":"number of features","group":label,"value":value,"statistic":"count","p_value":None,"p_adj":None,"n":n,"notes":"Source: Supplemental Tables 2-3 and manuscript summary."})
add("S19_Fig2","Table S19-2. Figure 2 source data","Compartment-level metabolite and protein counts used for Figure 2A.",pd.DataFrame(rows))

# Figure 3: metabolite categories and betaine validation.
rows=[]
for label,value in [("amino-acid derivatives",28),("anti-inflammatory mediators",12),("nucleotide derivatives",15)]:
    rows.append({"figure_panel":"Fig3A-C","x_variable":"metabolite class","y_variable":"candidate metabolite count","group":label,"value":value,"statistic":"count","p_value":None,"p_adj":None,"n":None,"notes":"Functional candidate categories."})
rows.append({"figure_panel":"Fig3D","x_variable":"treatment","y_variable":"total hemocyte count change","group":"betaine vs PBS","value":None,"statistic":"two-tailed unpaired Student t-test","p_value":0.028,"p_adj":None,"n":None,"notes":"Only betaine produced a significant increase among tested metabolites."})
add("S19_Fig3","Table S19-3. Figure 3 source data","Candidate metabolite category counts and betaine validation result.",pd.DataFrame(rows))

# Figure 4: DEP counts.
rows=[]
for tp,up,down in [("6 h",72,104),("12 h",73,101),("24 h",72,102)]:
    for d,v in [("up",up),("down",down)]: rows.append({"figure_panel":"Fig4A","x_variable":"time after hemospasia","y_variable":"number of DEPs","group":d,"value":v,"statistic":"count","p_value":None,"p_adj":None,"n":20,"notes":f"{tp} vs 0 h; source: Supplemental Table 5 and manuscript."})
rows.append({"figure_panel":"Fig4B","x_variable":"time-point intersection","y_variable":"core DEPs","group":"6/12/24 h shared","value":119,"statistic":"count","p_value":None,"p_adj":None,"n":20,"notes":"DEPs consistently regulated at all post-hemospasia time points."})
add("S19_Fig4","Table S19-4. Figure 4 source data","Differentially expressed protein counts and core temporal overlap.",pd.DataFrame(rows))

# Figure 5: single-cell atlas.
adata=ad.read_h5ad(A/"step3"/"step3_lyo_scores.h5ad"); obs=adata.obs.copy(); comp=pd.read_csv(A/"composition_new_names_10groups.csv"); rows=[]
for sample in ["blank","treated"]:
    d=obs[obs["sample"].astype(str).eq(sample)]
    rows.append({"figure_panel":"Fig5B","x_variable":"library","y_variable":"retained nuclei","group":SAMPLE_MAP[sample],"value":len(d),"statistic":"count","p_value":None,"p_adj":None,"n":len(d),"notes":f"Median genes={d.n_genes_by_counts.median():.0f}; median UMIs={d.total_counts.median():.0f}; median mito fraction={d.pct_counts_mt.median()/100:.4f}."})
for r in comp.itertuples(index=False):
    for grp,pct,n in [("control",r.pct_blank,r.n_blank),("hemospasia",r.pct_treated,r.n_treated)]:
        rows.append({"figure_panel":"Fig5G","x_variable":"state","y_variable":"composition (%)","group":f"{r.group}|{grp}","value":float(pct),"statistic":"library proportion","p_value":None,"p_adj":None,"n":int(n),"notes":"Composition source data; exploratory tests in Table S12."})
coords=adata.obsm["X_tsne"]
for i,cell in enumerate(adata.obs_names):
    grp=str(obs.iloc[i].cell_type); smp=SAMPLE_MAP[str(obs.iloc[i]["sample"])]
    rows.append({"figure_panel":"Fig5C,E,F","x_variable":"tSNE1","y_variable":"tSNE2","group":grp,"value":float(coords[i,0]),"statistic":"x coordinate","p_value":None,"p_adj":None,"n":1,"notes":f"nucleus_id={cell}; sample={smp}; x={coords[i,0]:.6f}; y={coords[i,1]:.6f}"})
    rows.append({"figure_panel":"Fig5C,E,F","x_variable":"tSNE1","y_variable":"tSNE2","group":grp,"value":float(coords[i,1]),"statistic":"y coordinate","p_value":None,"p_adj":None,"n":1,"notes":f"nucleus_id={cell}; sample={smp}; x={coords[i,0]:.6f}; y={coords[i,1]:.6f}"})
b=pd.read_csv(A/"Fig5_marker_bubble_10groups_data.csv")
for r in b.itertuples(index=False):
    rows.append({"figure_panel":"Fig5D","x_variable":"state","y_variable":"gene","group":str(r.group),"value":float(r.mean_logexpr),"statistic":"mean log-normalized expression","p_value":None,"p_adj":None,"n":None,"notes":f"gene={r.gene}; GeneID={r.gid}; pct_expressing={r.pct_expressing:.2f}; zscore={r.zscore:.3f}"})
add("S19_Fig5","Table S19-5. Figure 5 source data","Source data for the Hs atlas, QC, t-SNE, marker heat map and composition panels.",pd.DataFrame(rows))

# Figure 6: PH state source data.
rows=[]; cc=pd.read_csv(A/"figures_pub"/"2026-10-21_Fig2"/"preview"/"Fig2a_all_states_cellcycle.csv")
for r in cc.itertuples(index=False): rows.append({"figure_panel":"Fig6A","x_variable":"state","y_variable":"G2M (%)","group":str(r.cell_type),"value":float(r.G2M),"statistic":"cell-cycle phase fraction","p_value":None,"p_adj":None,"n":int(r.n),"notes":f"G1={r.G1:.2f}%; S={r.S:.2f}%; G2M={r.G2M:.2f}%"})
ann=pd.read_csv(A/"Fig5_dynamic_subcluster_annotation_table.csv"); ann=ann[(ann.parent=="PH")&(ann.key.isin(["PH_s0","PH_s1","PH_s2","PH_s3"]))]
for r in ann.itertuples(index=False):
    for grp,pct,n in [("control",r.pct_control,r.n_control),("hemospasia",r.pct_hemospasia,r.n_hemospasia)]: rows.append({"figure_panel":"Fig6B","x_variable":"PH substate","y_variable":"composition (%)","group":f"{r.key}|{grp}","value":float(pct),"statistic":"within-PH proportion","p_value":None,"p_adj":None,"n":int(n),"notes":f"label={r.proposed_state_label}"})
pc=pd.read_csv(A/"figures_pub"/"2026-10-23_Fig2D_PH_tSNE_v2"/"PH_only_tsne_coordinates.csv")
for r in pc.itertuples(index=False):
    cell=str(getattr(r,"_0",getattr(r,"Unnamed: 0","PH"))); grp=f"{r.group}|PH_s{int(r.subcluster)}"
    rows.append({"figure_panel":"Fig6D","x_variable":"tSNE1","y_variable":"tSNE2","group":grp,"value":float(r.tsne1),"statistic":"x coordinate","p_value":None,"p_adj":None,"n":1,"notes":f"nucleus_id={cell}; x={r.tsne1:.6f}; y={r.tsne2:.6f}"})
    rows.append({"figure_panel":"Fig6D","x_variable":"tSNE1","y_variable":"tSNE2","group":grp,"value":float(r.tsne2),"statistic":"y coordinate","p_value":None,"p_adj":None,"n":1,"notes":f"nucleus_id={cell}; x={r.tsne1:.6f}; y={r.tsne2:.6f}"})
uc=pd.read_csv(A/"step3"/"step3_ucell_by_celltype.csv")
for r in uc.itertuples(index=False):
    for prog,col in [("LyO_Progenitor","mean_LyO_Progenitor"),("LyO_Lymphocyte_precursor","mean_LyO_Lymphocyte_precursor")]: rows.append({"figure_panel":"Fig6E","x_variable":"state","y_variable":"mean UCell score","group":f"{r.cell_type}|{prog}","value":float(getattr(r,col)),"statistic":"UCell-style score","p_value":None,"p_adj":None,"n":int(r.n),"notes":"Mean score across the state."})
co=pd.read_csv(A/"cross_organ_progenitor"/"cross_organ_adjusted_summary.csv")
for r in co.itertuples(index=False): rows.append({"figure_panel":"Fig6F","x_variable":"dataset/target","y_variable":"adjusted progenitor potential (z)","group":str(r.target),"value":float(r.mean_adjusted_prog_z),"statistic":"cycle-stratified permutation","p_value":float(r.permutation_p),"p_adj":None,"n":int(r.n_target),"notes":f"95% CI {r.adjusted_ci_low:.3f} to {r.adjusted_ci_high:.3f}"})
add("S19_Fig6","Table S19-6. Figure 6 source data","Source data for PH cell-cycle composition, substate composition, PH t-SNE, UCell programs and cross-organ adjusted identity.",pd.DataFrame(rows))

(OUT/"tables_S19.json").write_text(json.dumps({"tables":TABLES},ensure_ascii=False,allow_nan=False),encoding="utf-8")
print("wrote",len(TABLES))
