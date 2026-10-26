p = r"D:\桌面\single cell\analysis_output\scripts\127_syncytium_tests.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace('_sym = raw.var["sym"].astype(str).values',
                  '_sym = (raw.var["symbol"] if "symbol" in raw.var.columns else raw.var["sym"]).astype(str).values')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
