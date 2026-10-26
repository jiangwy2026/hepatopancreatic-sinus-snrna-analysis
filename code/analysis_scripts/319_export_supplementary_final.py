# -*- coding: utf-8 -*-
"""Export final supplementary single panels and assembled figures."""

from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
import shutil
ROOT=Path(r'D:\桌面\single cell')
SUP=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary'
PREV=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_previews'
FINAL=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_final'
FINAL.mkdir(parents=True,exist_ok=True)
font=r'C:\Windows\Fonts\timesbd.ttf'

def load(p):
 im=Image.open(p).convert('RGBA'); bg=Image.new('RGBA',im.size,'white'); bg.alpha_composite(im); return bg.convert('RGB')
def contain(im,w,h):
 im=im.copy(); im.thumbnail((w,h),Image.Resampling.LANCZOS); return im
def save_raster(im,stem):
 im.save(stem.with_suffix('.png'),dpi=(300,300)); im.save(stem.with_suffix('.tiff'),dpi=(1200,1200),compression='tiff_lzw'); im.save(stem.with_suffix('.pdf'),'PDF',resolution=300,dpi=(300,300))

# Rebuild the wide LyO mirror as a square page from its three source panels.
orig=Image.open(SUP/'Supp_LyO_marker_mirror.png').convert('RGB'); d0=ImageDraw.Draw(orig)
for box in [(145,5,225,65),(1540,5,1630,65),(160,740,250,825)]: d0.rectangle(box,fill='white')
A=orig.crop((20,0,1410,720)); B=orig.crop((1510,0,2355,720)); C=orig.crop((140,750,2355,1490))
A2=contain(A,1450,760); B2=contain(B,720,620); C2=contain(C,1980,700)
SQ=2200; canvas=Image.new('RGB',(SQ,SQ),'white'); d=ImageDraw.Draw(canvas)
d.text((SQ//2,35),'Lymphoid-organ progenitor program mirror',font=ImageFont.truetype(font,38),fill='#0f172a',anchor='mm')
y=95; canvas.paste(A2,((SQ-A2.width)//2,y)); y+=A2.height+22; canvas.paste(B2,((SQ-B2.width)//2,y-(B2.height-520 if B2.height>520 else 0))); y+=B2.height+22; canvas.paste(C2,((SQ-C2.width)//2,max(y,SQ-50-C2.height)))
save_raster(canvas,FINAL/'S1_C_LyO_marker_mirror_square')

# Copy/re-export S1 A-B panels.
shutil.copy2(SUP/'Supp_PH_substate_markers_compact_legend.png',FINAL/'S1_A_PH_substate_markers.png')
shutil.copy2(SUP/'Supp_PH_substate_markers_compact_legend.svg',FINAL/'S1_A_PH_substate_markers.svg')
shutil.copy2(SUP/'Supp_PH_substate_markers_compact_legend.tiff',FINAL/'S1_A_PH_substate_markers.tiff')
shutil.copy2(SUP/'Supp_PH_substate_markers_compact_legend.pdf',FINAL/'S1_A_PH_substate_markers.pdf')
shutil.copy2(SUP/'Supp_PH_identity_vs_cellcycle.png',FINAL/'S1_B_PH_identity_vs_cellcycle.png')
shutil.copy2(SUP/'Supp_PH_identity_vs_cellcycle.svg',FINAL/'S1_B_PH_identity_vs_cellcycle.svg')
# TIFF/PDF for S1 B from PNG
im=Image.open(SUP/'Supp_PH_identity_vs_cellcycle.png').convert('RGB'); save_raster(im,FINAL/'S1_B_PH_identity_vs_cellcycle')

# S1 composite: A/B top, square C bottom, no added panel letters.
A=load(FINAL/'S1_A_PH_substate_markers.png'); B=load(FINAL/'S1_B_PH_identity_vs_cellcycle.png'); C=load(FINAL/'S1_C_LyO_marker_mirror_square.png')
W=2400; margin=40; gap=30; header=100; top_h=1000; colw=(W-2*margin-gap)//2
A2=contain(A,colw,top_h); B2=contain(B,colw,top_h); C2=contain(C,1850,1700)
y1=margin+header; y2=y1+top_h+gap; H=y2+C2.height+margin
S1=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(S1); d.text((W//2,35),'Supplementary Figure S1 | PH substate identity evidence',font=ImageFont.truetype(font,38),fill='#0f172a',anchor='mm')
S1.paste(A2,(margin+(colw-A2.width)//2,y1)); S1.paste(B2,(margin+colw+gap+(colw-B2.width)//2,y1)); S1.paste(C2,((W-C2.width)//2,y2))
save_raster(S1,FINAL/'Supplementary_Figure_S1')

# S2 individual panels: copy source formats and add PDF.
mapping={
 'S2_A_niche_marker_expression':('Supp_Niche_marker_dotplot.png','Supp_Niche_marker_dotplot.svg','Supp_Niche_marker_dotplot.tiff'),
 'S2_B_ligand_receptor_pairs':('Supp_Niche_LR_pairs.png','Supp_Niche_LR_pairs.svg','Supp_Niche_LR_pairs.tiff'),
 'S2_C_orthogroup_presence':('Supp_Niche_orthogroup_presence.png','Supp_Niche_orthogroup_presence.svg','Supp_Niche_orthogroup_presence.tiff'),
 'S2_D_receptor_expression':('Supp_Niche_receptor_expression.png','Supp_Niche_receptor_expression.svg','Supp_Niche_receptor_expression.tiff'),
 'S2_E_communication_sensitivity':('Supp_Niche_communication_sensitivity.png','Supp_Niche_communication_sensitivity.svg','Supp_Niche_communication_sensitivity.tiff')}
for stem,(png,svg,tif) in mapping.items():
 for src,ext in [(png,'.png'),(svg,'.svg'),(tif,'.tiff')]: shutil.copy2(SUP/src,FINAL/(stem+ext))
 im=Image.open(SUP/png).convert('RGB'); im.save(FINAL/(stem+'.pdf'),'PDF',resolution=300,dpi=(300,300))

# S2 readable composite from the earlier readability-optimized preview.
S2=load(PREV/'Supplementary_Figure_S2_niche_preview_readable.png'); save_raster(S2,FINAL/'Supplementary_Figure_S2_niche')
print(FINAL)
