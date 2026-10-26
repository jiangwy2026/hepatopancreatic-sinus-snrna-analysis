# -*- coding: utf-8 -*-
"""Resolve the orientation of paga transitions_confidence by correlating with a direct
group-level net-flow computed from the velocity graph (rows = source cells, verified)."""
import os, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, scvelo as scv
import matplotlib; matplotlib.use("Agg")
OUT = r"D:\桌面\single cell\analysis_output"
ORDER = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
A = ad.read_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"))
ct = A.obs["ct"].astype(str).values
VG = A.uns["velocity_graph"]; VGn = A.uns.get("velocity_graph_neg")
T = sp.csr_matrix(VG - VGn) if VGn is not None else sp.csr_matrix(VG)
T.data = np.clip(T.data, 0, None)
rs = np.asarray(T.sum(axis=1)).ravel(); rs[rs <= 0] = 1.0
P = sp.diags(1.0 / rs) @ T
gi = {g: i for i, g in enumerate(ORDER)}
src = np.array([gi[g] for g in ct]); ncell = np.bincount(src, minlength=9)
Tr = sp.csr_matrix((np.ones(len(ct)), (src, np.arange(len(ct)))), shape=(9, len(ct)))
flow = np.asarray((Tr @ P @ Tr.T).todense()) / ncell[:, None]
net_mine = flow - flow.T                       # rows = source (verified convention)
A.obs["ct"] = pd.Categorical(ct, categories=ORDER)
scv.tl.paga(A, groups="ct")
stored = A.uns["paga"]["transitions_confidence"]
stored = np.asarray(stored.todense()) if hasattr(stored, "todense") else np.asarray(stored)
off = ~np.eye(9, dtype=bool)
print("corr(net_mine, stored)     = %+.3f" % np.corrcoef(net_mine[off], stored[off])[0, 1])
print("corr(net_mine, stored.T)   = %+.3f" % np.corrcoef(net_mine[off], stored.T[off])[0, 1])
print("\nmy net flow (row -> col, positive = net outgoing), rounded:")
print(pd.DataFrame(np.where(net_mine > 0, net_mine, 0), index=ORDER, columns=ORDER).round(4).to_string())
print("\nstored transitions_confidence as printed before:")
print(pd.DataFrame(np.where(stored > 0, stored, 0), index=ORDER, columns=ORDER).round(4).to_string())
