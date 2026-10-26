# -*- coding: utf-8 -*-
"""91_velocity_flow.py - row-normalised group-to-group velocity flow + preview figures."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc, scvelo as scv

OUT = r"D:\桌面\single cell\analysis_output"
FIG = os.path.join(OUT, "figures_velocity"); os.makedirs(FIG, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORDER = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"))
    ct = A.obs["ct"].astype(str).values
    VG = A.uns["velocity_graph"]
    VGn = A.uns.get("velocity_graph_neg")
    T = sp.csr_matrix(VG - VGn) if VGn is not None else sp.csr_matrix(VG)
    T.data = np.clip(T.data, 0, None)
    rs = np.asarray(T.sum(axis=1)).ravel(); rs[rs <= 0] = 1.0
    P = sp.diags(1.0 / rs) @ T                      # row-normalised transition probabilities
    gi = {g: i for i, g in enumerate(ORDER)}
    src = np.array([gi[g] for g in ct])
    ncell = np.bincount(src, minlength=len(ORDER))
    Tr = sp.csr_matrix((np.ones(len(ct)), (src, np.arange(len(ct)))), shape=(len(ORDER), len(ct)))
    flow = np.asarray((Tr @ P @ Tr.T).todense())
    flow = flow / ncell[:, None]                    # mean outgoing probability per source cell
    fdf = pd.DataFrame(flow, index=ORDER, columns=ORDER)
    print("--- normalised group-to-group transition probability (row = source) ---")
    print(fdf.round(4).to_string())
    print("\nrow sums (should be ~1):", fdf.sum(axis=1).round(3).to_dict())
    iph = gi["PH"]
    print("\n--- directionality of the hemocyte lineage (PH as candidate root) ---")
    out = {}
    for g in ORDER:
        if g == "PH": continue
        ig = gi[g]
        f_fwd, f_rev = flow[iph, ig], flow[ig, iph]
        tot_i = flow[ig].sum()
        out[g] = {"PH_to_%s" % g: float(f_fwd), "%s_to_PH" % g: float(f_rev),
                  "net_PH_out": float(f_fwd - f_rev),
                  "PH_share_of_%s_outflow" % g: float(f_rev / tot_i) if tot_i else float("nan")}
        print("   PH->%-4s %.4f  (%.1f%% of PH outflow) | %-4s->PH %.4f (%.1f%% of %s outflow) | net(PH out) %+.4f" % (
            g, f_fwd, 100 * f_fwd / flow[iph].sum(), g, f_rev, 100 * f_rev / tot_i, g, f_fwd - f_rev))
    fdf.to_csv(os.path.join(OUT, "velocity_group_flow_normalized.csv"))
    json.dump(out, open(os.path.join(OUT, "VELOCITY_DIRECTIONALITY.json"), "w"), indent=1)

    scv.settings.figdir = FIG; sc.settings.figdir = FIG
    pal = [COL[g] for g in ORDER]
    scv.pl.velocity_embedding_stream(A, basis="umap", color="ct", palette=pal, legend_loc="on data",
                                     title="RNA velocity (stochastic), 9 functional groups",
                                     dpi=200, show=False, save="_stream.png")
    scv.pl.velocity_embedding(A, basis="umap", color="ct", palette=pal, dpi=200,
                              show=False, save="_arrows.png")
    sc.pl.umap(A, color=["velocity_pseudotime", "velocity_confidence"], cmap="viridis",
               save="_pseudotime.png", show=False)
    print("\nfigures ->", FIG)
    print("files:", sorted(os.listdir(FIG)))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\velocity_flow.log", "w", encoding="utf8", buffering=1)
    main()
