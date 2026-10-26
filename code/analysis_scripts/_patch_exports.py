import os, re
S = r"D:\桌面\single cell\analysis_output\scripts"

def read(p):
    return open(p, encoding="utf8").read().replace("\r\n", "\n")

def write(p, s):
    open(p, "w", encoding="utf8", newline="\n").write(s)

LOOP = '''{ind}bp = os.path.join(FIG, "{base}")
{ind}for _e, _kw in (("_preview.png", dict(dpi={dpi})), (".pdf", dict({extra})), (".svg", dict({extra})),
{ind}               (".tiff", dict({extra}dpi=1200, pil_kwargs={{"compression": "tiff_lzw"}}))):
{ind}    fig.savefig(bp + _e, **{ {kw} }
{ind}    print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
{ind}plt.close(fig)'''

def patch_generic(fn, base_new):
    p = os.path.join(S, fn); src = read(p)
    pat = re.compile(r'([ \t]*)p = os\.path\.join\(FIG, "([^"]+)_preview\.png"\)\n[ \t]*fig\.savefig\(p, ([^)]*)\); plt\.close\(fig\)')
    m = pat.search(src)
    assert m, fn
    ind, args = m.group(1), m.group(3)
    dpi = re.search(r"dpi=(\d+)", args)
    dpi = dpi.group(1) if dpi else "400"
    extra = ", ".join(a.strip() for a in args.split(",") if a.strip() and not a.strip().startswith("dpi="))
    new = (f'{ind}bp = os.path.join(FIG, "{base_new}")\n'
           f'{ind}for _e, _kw in (("_preview.png", dict(dpi={dpi})), (".pdf", dict({extra})), (".svg", dict({extra})),\n'
           f'{ind}               (".tiff", dict({extra}dpi=1200, pil_kwargs={{"compression": "tiff_lzw"}}))):\n'
           f'{ind}    fig.savefig(bp + _e, **{{**_kw}})\n'
           f'{ind}    print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))\n'
           f'{ind}plt.close(fig)')
    src = src[:m.start()] + new + src[m.end():]
    if extra:
        src = src.replace(f'{ind}    fig.savefig(bp + _e, **{{**_kw}})', f'{ind}    fig.savefig(bp + _e, {extra}, **_kw)')
    write(p, src); print("patched", fn, "| args:", args, "-> extra:", extra)

patch_generic("100_traj_lines_final.py", "Fig4E_lineages")
patch_generic("112_check_and_rerender_heatmap.py", "Fig2_marker_heatmap_singlecell_logcells")
patch_generic("119_figE_tsne_legend.py", "Fig4E_lineages_tsne")
patch_generic("120_tsne_gene_v4.py", "Fig_tsne_genefeatures")

# Fig3: single preview save -> four formats (base already named bp)
p = os.path.join(S, "73_fig3_final.py"); src = read(p)
old = 'fig.savefig(bp + "_preview.png", dpi=300, bbox_inches=None)\nplt.close()'
new = ('for _e, _kw in (("_preview.png", dict(dpi=300)), (".pdf", {}), (".svg", {}),\n'
       '               (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):\n'
       '    fig.savefig(bp + _e, bbox_inches=None, **_kw)\n'
       '    print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))\n'
       'plt.close()')
assert old in src, "73"
write(p, src.replace(old, new, 1)); print("patched 73_fig3_final.py")
