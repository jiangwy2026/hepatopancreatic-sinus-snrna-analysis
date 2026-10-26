p = r"D:\桌面\single cell\analysis_output\scripts\146_score_published_hep_markers.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace('sym = A.var["symbol"].astype(str).values; ids = A.var_names.astype(str).values',
                  'sym = (A.var["sym"] if "sym" in A.var.columns else A.var["symbol"]).astype(str).values; ids = A.var_names.astype(str).values')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
