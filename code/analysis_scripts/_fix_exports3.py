import os
S = r"D:\桌面\single cell\analysis_output\scripts"
for fn in ("112_check_and_rerender_heatmap.py", "120_tsne_gene_v4.py"):
    p = os.path.join(S, fn)
    src = open(p, encoding="utf8").read().replace("\r\n", "\n")
    assert 'print("saved", p)' in src, fn
    open(p, "w", encoding="utf8", newline="\n").write(src.replace('print("saved", p)', 'print("saved", bp + "_preview.png")'))
    print("fixed", fn)
