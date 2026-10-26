# -*- coding: utf-8 -*-
"""
89_velocity_scvelo.py - RNA velocity (scVelo) on the custom spliced/unspliced counts.
Question: does the velocity field run PH -> GH/MH1/MH2, and what are the per-group pseudotimes?
Windows note: multiprocessing (used by scv.tl.velocity_graph) needs the __main__ guard.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")          # must be set before pyplot/scvelo import (no Tk in this environment)
import numpy as np, pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scvelo as scv

ROOT = r"D:\桌面\single cell"; OUT = os.path.join(ROOT, "analysis_output")
FIG = os.path.join(OUT, "figures_velocity"); os.makedirs(FIG, exist_ok=True)
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B",
        "8": "NE", "9": "Neu"}
ORDER = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}

def main():
    scv.settings.verbosity = 2; scv.settings.n_jobs = 8
    scv.settings.figdir = FIG; sc.settings.figdir = FIG
    A = ad.read_h5ad(os.path.join(OUT, "adata_velocity_raw.h5ad"))
    H = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg.h5ad"))
    assert (A.obs_names == H.obs_names).all()
    ct = H.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    A.obs["ct"] = pd.Categorical(ct, categories=ORDER)
    A.obsm["X_umap"] = H.obsm["X_umap"].copy()
    print("groups:", pd.Series(ct).value_counts().reindex(ORDER).to_dict(), flush=True)

    scv.pp.filter_and_normalize(A, min_shared_counts=20)
    print("genes kept:", A.n_vars, flush=True)
    for _k in (None, "spliced", "unspliced", "ambiguous"):
        _M = A.X if _k is None else A.layers.get(_k)
        if _M is None: continue
        if sp.issparse(_M):
            _bad = ~np.isfinite(_M.data)
            if _bad.any():
                _M.data[_bad] = 0.0; print("  cleaned %d non-finite in %s" % (int(_bad.sum()), _k), flush=True)
        else:
            _bad = ~np.isfinite(_M)
            if _bad.any():
                _M[_bad] = 0.0; print("  cleaned %d non-finite in %s" % (int(_bad.sum()), _k), flush=True)
    scv.pp.moments(A, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(A, mode="stochastic")
    A.uns.pop("velocity_graph", None)
    scv.tl.velocity_graph(A, show_progress_bar=False)
    scv.tl.velocity_embedding(A, basis="umap")
    scv.tl.velocity_confidence(A)
    A.obs["is_root"] = pd.Series((ct == "PH").astype(bool), index=A.obs_names)
    print("is_root dtype:", A.obs["is_root"].dtype, "| n root:", int(A.obs["is_root"].sum()), flush=True)
    # velocity pseudotime with PH as the root (scvelo's wrapper trips over a pandas-2 categorical
    # check on Windows; drive the underlying VPT class directly instead)
    try:
        from scvelo.tools.velocity_pseudotime import VPT
        vpt = VPT(A.copy(), n_dcs=10, allow_kendall_tau_shift=True)
        _T = A.uns["velocity_graph"] - A.uns["velocity_graph_neg"]
        vpt._connectivities = _T + _T.T
        vpt.compute_transitions()
        vpt.compute_eigen(n_comps=10)
        vpt.set_iroot("is_root")   # column name (bool array is not accepted by scvelo)
        vpt.compute_pseudotime(inverse=False)
        A.obs["velocity_pseudotime"] = np.asarray(vpt.pseudotime)
        print("velocity_pseudotime (root = PH cells) computed", flush=True)
    except Exception as e:
        print("velocity_pseudotime failed (%r)" % (e,), flush=True)
        A.obs["velocity_pseudotime"] = np.nan

    try:
        scv.tl.latent_time(A)
    except Exception as e:
        print("latent_time unavailable (%r) - stochastic mode has no dynamical fit; "
              "reporting velocity_pseudotime only" % (e,), flush=True)
        A.obs["latent_time"] = np.nan
    print("velocity computed", A.shape, flush=True)

    lt = A.obs["latent_time"].values; vp = A.obs["velocity_pseudotime"].values
    conf = A.obs["velocity_confidence"].values
    rows = []
    for g in ORDER:
        m = ct == g
        rows.append({"group": g, "n": int(m.sum()), "latent_time": float(np.nanmean(lt[m])),
                     "velocity_pseudotime": float(np.nanmean(vp[m])),
                     "velocity_confidence": float(np.nanmean(conf[m]))})
    tab = pd.DataFrame(rows)
    print("\n--- per-group pseudotime / confidence ---", flush=True)
    print(tab.to_string(index=False), flush=True)

    VG = A.uns["velocity_graph"]
    VG = sp.csr_matrix(VG) if not sp.issparse(VG) else VG.tocsr()
    gi = {g: i for i, g in enumerate(ORDER)}
    src = np.array([gi[g] for g in ct])
    ncell_g = np.bincount(src, minlength=len(ORDER))
    Tr = sp.csr_matrix((np.ones(len(ct)), (src, np.arange(len(ct)))), shape=(len(ORDER), len(ct)))
    flow = np.asarray((Tr @ VG @ Tr.T).todense()) / np.maximum(ncell_g[:, None], 1)
    fdf = pd.DataFrame(flow, index=ORDER, columns=ORDER)
    print("\n--- mean velocity transition probability (row=source, col=target) ---", flush=True)
    print(fdf.round(4).to_string(), flush=True)

    iph = gi["PH"]; net = {}
    print("\n--- hemocyte-lineage net flow vs PH ---", flush=True)
    for g in ORDER:
        if g == "PH": continue
        ig = gi[g]
        d = flow[iph, ig] - flow[ig, iph]
        net[g] = {"PH_to": float(flow[iph, ig]), "to_PH": float(flow[ig, iph]), "net_PH_out": float(d)}
        print("   PH->%-4s %.4f | %-4s->PH %.4f | net(PH out) %+.4f" % (g, flow[iph, ig], g, flow[ig, iph], d), flush=True)

    A.write_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"))
    tab.to_csv(os.path.join(OUT, "velocity_per_group.csv"), index=False)
    fdf.to_csv(os.path.join(OUT, "velocity_group_flow.csv"))
    json.dump({"per_group": rows, "net_flow_vs_PH": net, "order": ORDER},
              open(os.path.join(OUT, "VELOCITY_SCVELO_REPORT.json"), "w"), indent=1)
    print("\nsaved adata_velocity_processed.h5ad", flush=True)

    scv.pl.velocity_embedding_stream(A, basis="umap", color="ct", palette=COL, legend_loc="on data",
                                     title="RNA velocity (stochastic) - 9 functional groups",
                                     save="_velocity_stream_preview.png", dpi=220, show=False)
    sc.pl.umap(A, color=["latent_time", "velocity_pseudotime", "velocity_confidence"],
               save="_velocity_pt_preview.png", show=False)
    scv.pl.velocity_embedding(A, basis="umap", color="ct", palette=COL,
                              save="_velocity_arrows_preview.png", dpi=220, show=False)
    print("figures in", FIG, flush=True)

if __name__ == "__main__":
    sys.stdout = open(os.path.join(ROOT, "logs", "velocity_scvelo.log"), "w", encoding="utf8", buffering=1)
    main()
