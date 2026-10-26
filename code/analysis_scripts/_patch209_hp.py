p=r"D:\桌面\single cell\analysis_output\scripts\209_second_fig_three_validations.py"
s=open(p,encoding="utf-8").read()
s=s.replace('def make_group_z(A, group_col, sym_map):','def make_group_z(A, group_col, sym_map, gene_names=None):')
s=s.replace('names = A.var["sym"].astype(str).tolist() if "sym" in A.var.columns else A.var_names.astype(str).tolist()','names = gene_names if gene_names is not None else (A.var["sym"].astype(str).tolist() if "sym" in A.var.columns else A.var_names.astype(str).tolist())')
s=s.replace('gg,gz,gnames,gmu,gsd=make_group_z(G,"cl",gsym)','hp_names=[hs_gid2sym.get(("ncbi_"+str(x)) if not str(x).startswith("ncbi_") else str(x), str(x)) for x in G.var_names]\n    gg,gz,gnames,gmu,gsd=make_group_z(G,"cl",gsym,gene_names=hp_names)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
