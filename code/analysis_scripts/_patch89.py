p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = """scv.settings.verbosity = 2
scv.settings.n_jobs = 8
"""
new = """scv.settings.verbosity = 2
scv.settings.n_jobs = 8
FIG = os.path.join(OUT, "figures_velocity")
os.makedirs(FIG, exist_ok=True)
scv.settings.figdir = FIG
sc.settings.figdir = FIG
scv.settings.set_figure_params("scvelo", dpi=200, fontsize=11)
"""
assert old in src
src = src.replace(old, new, 1)
src = src.replace('''scv.pl.velocity_embedding_stream(a, basis="X_umap", color="group", legend_loc="on data",
                                 save="_velocity_stream_preview.png", dpi=200, show=False)''',
'''scv.pl.velocity_embedding_stream(a, basis="X_umap", color="group", legend_loc="on data",
                                 title="RNA velocity (stochastic) - 9 functional groups",
                                 save="_velocity_stream_preview.png", dpi=220, show=False)''')
open(p, "w", encoding="utf8").write(src)
print("patched 89")
