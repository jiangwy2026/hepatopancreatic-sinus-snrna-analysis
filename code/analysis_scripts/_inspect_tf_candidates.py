import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
R=pd.read_csv(os.path.join(O,"PH_TF_candidates.csv"))
T=pd.read_csv(os.path.join(O,"PH_TF_target_programs.csv"))
E=pd.read_csv(os.path.join(O,"PH_TF_target_enrichment_up.csv"))
for tf in ["Kdm4B","Cfp1","kn","br","pan","Taf7","foxo","Trl","Stat5b","pnt","lola","ct","osa","Adar","Tbp"]:
    r=R[R.TF==tf]
    e=E[E.TF==tf]
    print("\n",tf)
    if len(r): print(r[["TF","TF_family","early_logFC","early_pvals_adj","rho_G2M","rho_ccnd2","rho_pseudotime","target_rho_G2M","target_rho_ccnd2","target_diff_early"]].round(3).to_string(index=False))
    if len(e): print(e[["TF","overlap_up","hypergeom_p"]].round(3).to_string(index=False))
