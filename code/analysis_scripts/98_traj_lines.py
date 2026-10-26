# -*- coding: utf-8 -*-
"""98_traj_lines.py - add Slingshot-style lineage curves on the UMAP (diagnostics first)."""
import os, sys, json
import numpy as np, pandas as pd
import anndata as ad
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.interpolate import UnivariateSpline

OUT = r"D:\桌面\single cell\analysis_output"
ORD = ["PH", "GH", "MH1", "MH2"]

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"]); pt = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    cent = np.vstack([U[ct == g].mean(0) for g in ORD])
    print("centroids (UMAP):")
    for g, c in zip(ORD, cent):
        print("   %-4s %s  n=%d  meanPT=%.3f" % (g, np.round(c, 2).tolist(), int((ct == g).sum()), pt[ct == g].mean()))

    # MST over cluster centroids
    d = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2)
    np.fill_diagonal(d, 0)
    mst = minimum_spanning_tree(d).toarray()
    edges = [(ORD[i], ORD[j], round(d[i, j], 2)) for i in range(4) for j in range(4) if mst[i, j] > 0]
    print("\nMST edges over centroids:", edges)

    # adjacency, root = PH, lineages = paths from PH to every other node
    adj = {g: [] for g in ORD}
    for a, b, _ in edges:
        adj[a].append(b); adj[b].append(a)
    def paths(node, seen):
        out = []
        for nb in adj[node]:
            if nb in seen: continue
            p = paths(nb, seen | {nb})
            if not p: out.append([node, nb])
            else:
                for q in p: out.append([node] + q)
        return out
    lin = paths("PH", {"PH"})
    print("lineages from PH:", lin)
    json.dump({"centroids": {g: cent[i].tolist() for i, g in enumerate(ORD)},
               "mst_edges": edges, "lineages": lin},
              open(os.path.join(OUT, "TRAJ_LINEAGES.json"), "w"), indent=1)

    # fit a principal-curve-like smooth path for each lineage (bin mean positions along pseudotime)
    out = {}
    for path in lin:
        m = np.isin(ct, path)
        sub_pt = pt[m]; sub_U = U[m]
        o = np.argsort(sub_pt); sub_pt = sub_pt[o]; sub_U = sub_U[o]
        q = np.quantile(sub_pt, np.linspace(0, 1, 26))
        xs, ys, w = [], [], []
        for k in range(25):
            sel = (sub_pt >= q[k]) & (sub_pt <= q[k + 1])
            if sel.sum() < 15: continue
            xs.append(sub_U[sel, 0].mean()); ys.append(sub_U[sel, 1].mean()); w.append(sel.sum())
        t = np.arange(len(xs), dtype=float)
        spl_x = UnivariateSpline(t, np.array(xs), k=3, s=len(t) * 0.05)
        spl_y = UnivariateSpline(t, np.array(ys), k=3, s=len(t) * 0.05)
        tt = np.linspace(0, len(t) - 1, 200)
        curve = np.c_[spl_x(tt), spl_y(tt)]
        name = "->".join(path)
        out[name] = curve.tolist()
        print("lineage %-14s bins=%d  curve start %s end %s" % (
            name, len(t), np.round(curve[0], 2).tolist(), np.round(curve[-1], 2).tolist()))
    json.dump(out, open(os.path.join(OUT, "TRAJ_CURVES.json"), "w"))
    print("\nsaved TRAJ_CURVES.json")

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_lines.log", "w", encoding="utf8", buffering=1)
    main()
