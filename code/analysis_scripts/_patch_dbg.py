p = r"D:\桌面\single cell\analysis_output\scripts\_dbg_auc.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace('A.var["sym"][j]', 'A.var["sym"].values[j]')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
