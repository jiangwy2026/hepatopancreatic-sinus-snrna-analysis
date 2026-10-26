import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
T=pd.read_csv(os.path.join(O,"PH_TF_target_programs.csv"))
print(T.sort_values("target_rho_G2M",ascending=False).head(30).round(3).to_string(index=False))
print("\nBy ccnd2")
print(T.sort_values("target_rho_ccnd2",ascending=False).head(30).round(3).to_string(index=False))
