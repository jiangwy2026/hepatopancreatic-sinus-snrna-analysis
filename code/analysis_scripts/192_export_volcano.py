# -*- coding: utf-8 -*-
"""192_export_volcano.py - export standalone early-PH volcano as SVG + 1200-dpi TIFF."""
import os, sys, importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = r"D:\桌面\single cell"
BASE = os.path.join(ROOT, "analysis_output")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-17", "机制_火山图")
os.makedirs(DELIV, exist_ok=True)
spec = importlib.util.spec_from_file_location("m188", os.path.join(BASE, "scripts", "188_ccnd2_TGM_mechanism_preview.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

fig = plt.figure(figsize=(4.35, 4.05), facecolor="white")
ax = fig.add_axes([0.13, 0.14, 0.82, 0.79])
m.draw_volcano(ax)
stem = os.path.join(DELIV, "Fig_earlyPH_volcano")
svg = stem + ".svg"; tif = stem + ".tiff"
fig.savefig(svg, bbox_inches="tight", pad_inches=0.04, facecolor="white")
fig.savefig(tif, dpi=1200, bbox_inches="tight", pad_inches=0.04, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
with Image.open(tif) as im:
    if "A" in im.getbands():
        rgba = im.convert("RGBA"); bg = Image.new("RGB", rgba.size, (255,255,255))
        bg.paste(rgba, mask=rgba.getchannel("A")); tmp = tif + ".tmp.tiff"
        bg.save(tmp, format="TIFF", dpi=(1200,1200), compression="tiff_lzw")
        os.replace(tmp, tif)
plt.close(fig)
print("saved", svg)
print("saved", tif)
