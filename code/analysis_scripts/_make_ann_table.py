import os, pandas as pd
OUT = r"D:\桌面\single cell\analysis_output"
CL = ["0","1","2","3","4","5","6","7","8","9"]
NAME = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NEU","9":"NEU"}
E = pd.read_csv(os.path.join(OUT,"elife_label_transfer_z.csv"), index_col=0); E.columns=[str(c) for c in E.columns]
G = pd.read_csv(os.path.join(OUT,"geo_label_transfer_z.csv"), index_col=0); G.columns=[str(c) for c in G.columns]
mk = pd.read_csv(os.path.join(OUT,"marker_heatmap_selected_fixed.csv"))
rows=[]
for c in CL:
    e=E[c]; g=G[c]; m=mk[mk.grp==NAME[c]]
    rows.append({"cluster":"C"+c,"our_name":NAME[c],"n_markers_tested":len(m),
                 "top_markers":", ".join(m["gene"].astype(str).tolist()[:5]),
                 "best_eLife":e.idxmax(),"z_eLife":round(e.max(),2),
                 "best_GEO":g.idxmax(),"z_GEO":round(g.max(),2)})
T=pd.DataFrame(rows); print(T.to_string(index=False)); T.to_csv(os.path.join(OUT,"annotation_evidence_table.csv"), index=False)
