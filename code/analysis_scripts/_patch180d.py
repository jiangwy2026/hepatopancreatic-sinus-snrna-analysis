p = r"D:\桌面\single cell\analysis_output\scripts\180_pseudotime_tsne_all_subtypes_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('PH-rooted pseudotime (all subtypes)', 'Pseudotime (all subtypes)')
s = s.replace('        lab = "PH (root)" if g == "PH" else g', '        lab = g')
open(p,"w",encoding="utf-8").write(s)
print("patched")
