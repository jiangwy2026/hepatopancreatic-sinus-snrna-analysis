p = r"D:\桌面\single cell\analysis_output\scripts\148_specificity_figure.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace("        tot_ = mean.mean(1)\n        mc = int(np.argmax(tot_))", "        mc = int(np.argmax(mean.mean(axis=0)))")
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
