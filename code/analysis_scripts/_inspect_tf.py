import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
R=pd.read_csv(os.path.join(O,"PH_TF_candidates.csv"))
print("rows",len(R),"with target",R.target_rho_G2M.notna().sum())
print("\nBy target G2M")
print(R.dropna(subset=["target_rho_G2M"]).sort_values("target_rho_G2M",ascending=False)[["TF","TF_family","early_logFC","rho_G2M","rho_ccnd2","rho_pseudotime","target_rho_G2M","target_rho_ccnd2","target_diff_early"]].head(25).round(3).to_string(index=False))
print("\nBy rho_G2M")
print(R.sort_values("rho_G2M",ascending=False)[["TF","TF_family","early_logFC","rho_G2M","rho_ccnd2","rho_pseudotime","target_rho_G2M","target_rho_ccnd2"]].head(25).round(3).to_string(index=False))
print("\nBy rho_ccnd2")
print(R.sort_values("rho_ccnd2",ascending=False)[["TF","TF_family","early_logFC","rho_G2M","rho_ccnd2","rho_pseudotime","target_rho_G2M","target_rho_ccnd2"]].head(25).round(3).to_string(index=False))
