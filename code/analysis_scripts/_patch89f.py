p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = """    try:
        scv.tl.velocity_pseudotime(A, vkey="velocity", root_key="is_root")
    except Exception as e:
        print("velocity_pseudotime failed (%r); falling back to latent_time only" % (e,), flush=True)
        A.obs["velocity_pseudotime"] = np.nan"""
new = """    # velocity pseudotime with PH as the root (scvelo's wrapper trips over a pandas-2 categorical
    # check on Windows; drive the underlying VPT class directly instead)
    try:
        from scvelo.tools.velocity_pseudotime import VPT
        vpt = VPT(A.copy(), n_dcs=10, allow_kendall_tau_shift=True)
        _T = A.uns["velocity_graph"] - A.uns["velocity_graph_neg"]
        vpt._connectivities = _T + _T.T
        vpt.compute_transitions()
        vpt.compute_eigen(n_comps=10)
        vpt.set_iroot(A.obs["is_root"].values.astype(bool))
        vpt.compute_pseudotime(inverse=False)
        A.obs["velocity_pseudotime"] = np.asarray(vpt.pseudotime)
        print("velocity_pseudotime (root = PH cells) computed", flush=True)
    except Exception as e:
        print("velocity_pseudotime failed (%r)" % (e,), flush=True)
        A.obs["velocity_pseudotime"] = np.nan
"""
assert old in src
src = src.replace(old, new, 1)
old2 = "    scv.tl.latent_time(A)"
new2 = """    try:
        scv.tl.latent_time(A)
    except Exception as e:
        print("latent_time unavailable (%r) - stochastic mode has no dynamical fit; "
              "reporting velocity_pseudotime only" % (e,), flush=True)
        A.obs["latent_time"] = np.nan"""
assert old2 in src
src = src.replace(old2, new2, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
