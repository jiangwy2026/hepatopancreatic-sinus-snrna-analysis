p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = """    A.obs["is_root"] = (ct == "PH")
    scv.tl.velocity_pseudotime(A, vkey="velocity", root_key="is_root")"""
new = """    A.obs["is_root"] = pd.Series((ct == "PH").astype(bool), index=A.obs_names)
    print("is_root dtype:", A.obs["is_root"].dtype, "| n root:", int(A.obs["is_root"].sum()), flush=True)
    try:
        scv.tl.velocity_pseudotime(A, vkey="velocity", root_key="is_root")
    except Exception as e:
        print("velocity_pseudotime failed (%r); falling back to latent_time only" % (e,), flush=True)
        A.obs["velocity_pseudotime"] = np.nan"""
assert old in src
src = src.replace(old, new, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
