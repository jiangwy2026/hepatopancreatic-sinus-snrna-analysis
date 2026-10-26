# -*- coding: utf-8 -*-
"""226_assemble_second_figure_preview.py
Assemble a 6-panel preview for the treatment-response figure."""
from PIL import Image, ImageOps, ImageDraw, ImageFont
import os
ROOT=r"D:\桌面\single cell"; PARTS=os.path.join(ROOT,'analysis_output','figures_pub','fig5_parts'); DELIV=os.path.join(ROOT,'文章配图','最终_2026-10-20','第二张图_抽血处理细胞变化')
os.makedirs(DELIV,exist_ok=True)
panels=[
 ('a','Fig5_dynamic_composition_cellcycle_preview.png','Composition and cell cycle'),
 ('b','Fig5_parent_program_effect_heatmap_preview.png','QC-adjusted program changes'),
 ('c','Fig5_qc_standardized_composition_preview.png','QC-standardized composition'),
 ('d','Fig5_R_PMlike_substate_focus_preview.png','R and PM-like substates'),
 ('e','Fig5_PH_DPT_pseudotime_density_preview.png','PH DPT distribution'),
 ('f','Fig5_PH_DPT_pseudotime_heatmap_preview.png','PH pseudotime heatmap'),
]
CW,CH=1800,1150; pad=34; header=130
canvas=Image.new('RGB',(CW*2,header+CH*3),(255,255,255))
draw=ImageDraw.Draw(canvas)
try:
 font=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',44)
 small=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',34)
except Exception:
 font=ImageFont.load_default(); small=font
title='Hemospasia-induced hemocyte state remodeling'
tw=draw.textlength(title,font=font); draw.text(((canvas.width-tw)/2,26),title,fill=(20,20,20),font=font)
for i,(label,fn,caption) in enumerate(panels):
 r=i//2; c=i%2; x0=c*CW; y0=header+r*CH
 im=Image.open(os.path.join(PARTS,fn)).convert('RGB')
 box=(CW-2*pad,CH-2*pad)
 im=ImageOps.contain(im,box,Image.Resampling.LANCZOS)
 px=x0+(CW-im.width)//2; py=y0+(CH-im.height)//2
 canvas.paste(im,(px,py))
 draw.text((x0+18,y0+12),label,fill=(0,0,0),font=font)
 draw.text((x0+74,y0+20),caption,fill=(70,70,70),font=small)
out_png=os.path.join(PARTS,'Fig5_second_figure_composite_preview.png'); canvas.save(out_png,format='PNG',dpi=(300,300))
out_tif=os.path.join(DELIV,'Fig5_second_figure_composite_preview.tiff'); canvas.save(out_tif,format='TIFF',dpi=(600,600),compression='tiff_lzw')
print('saved',out_png); print('saved',out_tif)
