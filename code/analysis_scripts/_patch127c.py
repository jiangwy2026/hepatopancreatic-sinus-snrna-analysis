p = r"D:\桌面\single cell\analysis_output\scripts\127_syncytium_tests.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace("norm = sp.diags(1.0 / np.maximum(tot, 1)) @ X", "norm = sp.diags(1e4 / np.maximum(tot, 1)) @ X")
src = src.replace("Dp.sample.values", 'Dp["sample"].values')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
