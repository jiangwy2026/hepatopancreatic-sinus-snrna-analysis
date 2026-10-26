import os
S = r"D:\桌面\single cell\analysis_output\scripts"
def patch(fn, old, new):
    p = os.path.join(S, fn); src = open(p, encoding="utf8").read().replace("\r\n", "\n")
    assert old in src, fn
    open(p, "w", encoding="utf8", newline="\n").write(src.replace(old, new, 1)); print("patched", fn)

# 102 = 单面板 E（UMAP 版）-> 四格式
patch("102_figE_only.py",
 'p = os.path.join(FIG, "Fig4E_lineages_preview.png")\n    fig.savefig(p, dpi=400, bbox_inches="tight"); plt.close(fig)',
 '''bp = os.path.join(FIG, "Fig4E_lineages")
    for _e, _kw in (("_preview.png", dict(dpi=400)), (".pdf", {}), (".svg", {}),
                    (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + _e, bbox_inches="tight", **_kw)
        print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
    plt.close(fig)''')

# 100 = D+E 双面板 -> 改成 Fig4DE_ 前缀，避免与单面板 E 重名
patch("100_traj_lines_final.py",
      'bp = os.path.join(FIG, "Fig4E_lineages")',
      'bp = os.path.join(FIG, "Fig4DE_traj_pseudotime")')
