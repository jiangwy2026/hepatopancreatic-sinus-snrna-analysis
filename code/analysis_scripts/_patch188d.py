p=r"D:\桌面\single cell\analysis_output\scripts\188_ccnd2_TGM_mechanism_preview.py"
s=open(p,encoding="utf-8").read()
s=s.replace('''        r = D[D.gene == gene]
        if r.empty: continue
        r = r.iloc[0]''','''        r = D[D.gene == gene]
        if r.empty: continue
        r = r.sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False]).iloc[0]''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
