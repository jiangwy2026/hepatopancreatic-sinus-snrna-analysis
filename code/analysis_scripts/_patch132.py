p = r"D:\桌面\single cell\analysis_output\scripts\132_pertype_qc.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace("target_size=1e4", "target_sum=1e4")
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
