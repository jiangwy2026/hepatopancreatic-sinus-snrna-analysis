p = r"D:\桌面\single cell\analysis_output\scripts\149_percell_labeltransfer.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace('pd.Series(cl, name="our")', 'pd.Series(cl, index=A.obs_names, name="our")')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched", src.count('index=A.obs_names, name="our"'))
