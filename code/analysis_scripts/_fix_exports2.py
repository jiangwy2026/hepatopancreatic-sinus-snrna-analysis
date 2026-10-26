import os, re
S = r"D:\桌面\single cell\analysis_output\scripts"
for fn in ("112_check_and_rerender_heatmap.py", "119_figE_tsne_legend.py", "120_tsne_gene_v4.py"):
    p = os.path.join(S, fn)
    src = open(p, encoding="utf8").read().replace("\r\n", "\n")
    new = re.sub(r'fig\.savefig\(bp \+ _e, bbox_inches="tight", pad_inches=[0-9.]+, \*\*_kw\)',
                 'fig.savefig(bp + _e, **_kw)', src)
    assert new != src, fn
    open(p, "w", encoding="utf8", newline="\n").write(new)
    print("fixed", fn)
