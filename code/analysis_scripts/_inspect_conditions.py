import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
C=pd.read_csv(os.path.join(O,"Fig6_condition_program_scores.csv"))
print(C.head(30).to_string(index=False))
print("\nDatasets",C.dataset.unique())
# show real condition rows and p rows separately
for ds in C.dataset.unique():
 print("\n##",ds)
 d=C[(C.dataset==ds)&(~C.condition.str.endswith('__p_kruskal'))]
 print(d.pivot(index="program",columns="condition",values="mean_score").round(2).to_string())
 print("p",C[(C.dataset==ds)&(C.condition.str.endswith('__p_kruskal'))][["program","condition","mean_score"]].round(4).to_string(index=False))
