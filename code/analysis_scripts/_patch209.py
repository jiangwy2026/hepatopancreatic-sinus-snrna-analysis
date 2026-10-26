p=r"D:\桌面\single cell\analysis_output\scripts\209_second_fig_three_validations.py"
s=open(p,encoding="utf-8").read()
# make group names symbols
s=s.replace('return gs, Z, A.var_names.astype(str).tolist(), mu, sd','names = A.var["sym"].astype(str).tolist() if "sym" in A.var.columns else A.var_names.astype(str).tolist()\n    return gs, Z, names, mu, sd')
# map_set accept dict
s=s.replace('def map_set(genes, names, hs_gid2sym):','def map_set(genes, name2idx, hs_gid2sym):')
s=s.replace('        if s in names: out.append(names[s])','        if s in name2idx: out.append(name2idx[s])')
# all calls update names->name2idx
s=s.replace('idx=map_set(marker[c],name2idx,gid2sym)','idx=map_set(marker[c],name2idx,gid2sym)')
# condition_scores signature and body
s=s.replace('def condition_scores(A, group_col, cond_map, label):','def condition_scores(A, group_col, cond_map, label, gene_names):')
s=s.replace('X=A.X.tocsc(); names=A.var_names.astype(str).tolist(); name2idx={s:i for i,s in enumerate(names)}','X=A.X.tocsc(); names=gene_names; name2idx={s:i for i,s in enumerate(names)}')
# condition calls: pass symbol names
s=s.replace('condition_scores(H,"cond",lambda x:x,"Hepatopancreatic sinus (Hs)")','condition_scores(H,"cond",lambda x:x,"Hepatopancreatic sinus (Hs)", H.var["sym"].astype(str).tolist())')
s=s.replace('condition_scores(E,"cond",lambda x:x,"Hemolymph")','condition_scores(E,"cond",lambda x:x,"Hemolymph", E.var["sym"].astype(str).tolist())')
s=s.replace('condition_scores(G,"cond",lambda x:hp_cond.get(x,x),"Hepatopancreas")','condition_scores(G,"cond",lambda x:hp_cond.get(x,x),"Hepatopancreas", [gid2sym.get("ncbi_"+str(x), gid2sym.get(str(x), str(x))) for x in G.var_names])')
open(p,"w",encoding="utf-8").write(s)
print("patched")
