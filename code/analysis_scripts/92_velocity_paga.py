# -*- coding: utf-8 -*-
"""92_velocity_paga.py - size-corrected directed flow between groups (scVelo PAGA transition confidence)."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd
import anndata as ad
import scanpy as sc, scvelo as scv

OUT = r"D:\桌面\single cell\analysis_output"
FIG = os.path.join(OUT, "figures_velocity")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORDER = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}

def main():
    scv.settings.figdir = FIG; sc.settings.figdir = FIG
    A = ad.read_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"))
    A.obs["ct"] = pd.Categorical(A.obs["ct"].astype(str), categories=ORDER)
    scv.tl.paga(A, groups="ct")
    tc = A.uns["paga"]["transitions_confidence"]
    tc = pd.DataFrame(np.asarray(tc.todense()) if hasattr(tc, "todense") else np.asarray(tc),
                      index=ORDER, columns=ORDER)
    print("--- PAGA velocity transition confidence (row = source cluster) ---")
    print(tc.round(4).to_string())
    iph = ORDER.index("PH")
    print("\n--- outgoing edges from PH (size-corrected) ---")
    for g in ORDER:
        if g == "PH": continue
        print("   PH->%-4s  %.4f" % (g, tc.iloc[iph][g]))
    print("\n--- top 12 directed edges overall ---")
    st = tc.stack().sort_values(ascending=False)
    for (a, b), v in st.head(12).items():
        if a != b:
            print("   %-4s -> %-4s %.4f  (reverse %.4f)" % (a, b, v, tc.iloc[ORDER.index(b)][a]))
    tc.to_csv(os.path.join(OUT, "velocity_paga_transitions.csv"))
    json.dump({"paga_transitions": tc.to_dict()}, open(os.path.join(OUT, "VELOCITY_PAGA.json"), "w"), indent=1)

    scv.pl.paga(A, basis="umap", size=0.3, alpha=0.1, min_edge_width=1.5, node_size_scale=2.0,
                title="PAGA velocity graph (9 functional groups)", show=False,
                save="_paga_velocity.png")
    print("\nfigures:", sorted(os.listdir(FIG)))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\velocity_paga.log", "w", encoding="utf8", buffering=1)
    main()
