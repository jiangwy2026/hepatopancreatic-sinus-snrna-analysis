p=r"D:\桌面\single cell\analysis_output\scripts\183_y_branch_diffmap_elife_preview.py"
s=open(p,encoding="utf-8").read()
if 'from PIL import Image' not in s:
    s=s.replace('import numpy as np, pandas as pd, scanpy as sc, anndata as ad','import numpy as np, pandas as pd, scanpy as sc, anndata as ad\nfrom PIL import Image')
old='''    p = os.path.join(PREV, "02d_branch_diffmap_elife_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)'''
new='''    p = os.path.join(PREV, "02d_branch_diffmap_elife_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    # final exports
    deliv = os.path.join(r"D:\\桌面\\single cell\\文章配图\\最终_2026-10-17", "拟时分叉图")
    os.makedirs(deliv, exist_ok=True)
    stem = os.path.join(deliv, "Fig_pseudotime_branch_diffmap")
    fig.savefig(stem + ".svg", bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(stem + ".tiff", dpi=1200, bbox_inches="tight", pad_inches=0.04,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(stem + ".tiff") as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA"); bg = Image.new("RGB", rgba.size, (255,255,255))
            bg.paste(rgba, mask=rgba.getchannel("A")); tmp = stem + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(1200,1200), compression="tiff_lzw")
            os.replace(tmp, stem + ".tiff")
    print("saved", stem + ".svg")
    print("saved", stem + ".tiff")
    plt.close(fig)'''
s=s.replace(old,new)
open(p,"w",encoding="utf-8").write(s)
print("patched")
