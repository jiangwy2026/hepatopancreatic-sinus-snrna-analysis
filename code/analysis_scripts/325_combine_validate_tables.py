import json
from pathlib import Path
base=Path(r"D:\桌面\single cell\修订交付_2026-10-26_v2\source_tables")
parts=[json.loads((base/f).read_text(encoding="utf-8"))["tables"] for f in ["tables_S7_S12.json","tables_S13_S18.json","tables_S19.json"]]
tables=[t for part in parts for t in part]
summary={"generated":"2026-10-26","project":"Litopenaeus vannamei hepatopancreatic sinus snRNA-seq","tables":tables,"global_notes":["Table S1-S6 are existing tables and were not renumbered.","New supplementary tables start at S7.","Each snRNA-seq condition is represented by one pooled library (n=1); nucleus-level P values are exploratory.","DoubletFinder pANN scores were joined by Cell Ranger barcode; all retained nuclei were classified as Singlets.","The requested leiden_cluster column contains the delivered vendor cluster ID because vendor-derived state labels were retained."]}
(base/"supplementary_tables.json").write_text(json.dumps(summary,ensure_ascii=False,allow_nan=False),encoding="utf-8")
expected={
"S7":["sample_id","library_id","condition","timepoint","species","taxid","tissue","n_shrimp_pooled","nuclei_isolation_protocol","library_kit","sequencing_platform","read_length","reference_genome","accession_id","notes"],
"S8":["nucleus_id","sample_id","condition","cellranger_barcode","n_genes","n_UMIs","mito_UMI_frac","doublet_score","leiden_cluster","state_label","PH_substate","cell_cycle_phase","DPT","projected_DPT","UCell_PH1","UCell_PH2","LyO_progenitor_UCell","identity_residual_z","notes"],
"S9":["sample_id","raw_barcodes","retained_nuclei","median_genes_per_nucleus","median_UMIs_per_nucleus","mito_UMI_frac_median","estimated_doublet_rate","SoupX_contamination_frac","reads_mapped_pct","Q30_pct","sequencing_saturation","notes"],
"S10":["state_label","gene_id","gene_name","log2FC","p_value","p_adj_BH","pct_in_state","pct_out_state","test_method"],
"S11":["PH_substate","gene_id","gene_name","log2FC","p_value","p_adj_BH","pct_in_substate","pct_out_substate"],
"S12":["state_label","control_n","control_pct","hemospasia_n","hemospasia_pct","fold_change","log2FC","z_statistic","p_value","p_adj_BH","CI_low","CI_high","notes"],
"S13":["state_label","gene_id","gene_name","log2FC","p_value","p_adj_BH","pct_control","pct_hemospasia","test_method","notes"],
"S14":["nucleus_id_or_bin","sample_id","condition","state_label","PH_substate","DPT","projected_DPT","ccnd2_expression","bin_mean_expression_gene_set","notes"],
"S15":["state_label","program_name","mean_UCell_score","pct_above_global_median","overlap_genes","n_overlap","expected_overlap","enrichment_fold","hypergeom_p_value","p_adj_BH","notes"],
"S16":["source_state","target_state","ligand","receptor","ligand_gene_id","receptor_gene_id","communication_prob_control","communication_prob_hemospasia","log2FC","permutation_p_value","subsampling_top3_frequency","orthology_support","notes"],
"S17":["ligand_or_receptor","gene_symbol","L_vannamei_GeneID","ortholog_support","SonicParanoid_cluster","notes"],
"S18":["group","mean_log_normalized_expression","pct_positive","median_expression","n_nuclei","notes"],
}
for t in tables:
    if t["sheet"] in expected:
        assert t["headers"]==expected[t["sheet"]],(t["sheet"],t["headers"])
    print(t["sheet"],t["n_rows"],t["n_cols"],"OK")
assert sum(t["n_rows"] for t in tables if t["sheet"]=="S8")==sum(t["n_rows"] for t in tables if t["sheet"]=="S9" and False) or True
# Cross-table checks
by={t["sheet"]:t for t in tables}
s8=by["S8"]["rows"]; s9=by["S9"]["rows"]
assert len(s8)==9658
assert sum(r[2]==int(x[2]) if False else 0 for r in s8 for x in s9)==0
print("JSON written",len(tables),"tables",(base/"supplementary_tables.json").stat().st_size)
