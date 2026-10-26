p = r"D:\桌面\single cell\analysis_output\scripts\_resolve_paga_orientation.py"
src = open(p, encoding="utf8").read()
src = src.replace('stored = A.uns["paga"]["transitions_confidence"]',
                  'A.obs["ct"] = pd.Categorical(ct, categories=ORDER)\nscv.tl.paga(A, groups="ct")\nstored = A.uns["paga"]["transitions_confidence"]')
open(p, "w", encoding="utf8").write(src)
print("patched")
