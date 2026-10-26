from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
ROOT=Path(r'D:\桌面\single cell'); SUP=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary'; OUT=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_previews'
def load(p):
 im=Image.open(p).convert('RGBA'); bg=Image.new('RGBA',im.size,'white'); bg.alpha_composite(im); return bg.convert('RGB')
def contain(im,w,h):
 im=im.copy(); im.thumbnail((w,h),Image.Resampling.LANCZOS); return im
A=load(SUP/'Supp_PH_substate_markers_compact_legend.png'); B=load(SUP/'Supp_PH_identity_vs_cellcycle.png'); C=load(SUP/'Supp_LyO_marker_mirror_nolabels_full.png')
W=2400; margin=40; gap=30; header=100; top_h=1000; colw=(W-2*margin-gap)//2
A2=contain(A,colw,top_h); B2=contain(B,colw,top_h); C2=contain(C,W-2*margin,1500)
y1=margin+header; y2=y1+top_h+gap; H=y2+C2.height+margin
canvas=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(canvas); font=ImageFont.truetype(r'C:\Windows\Fonts\timesbd.ttf',38)
d.text((W//2,35),'Supplementary Figure S1 | PH substate identity evidence',font=font,fill='#0f172a',anchor='mm')
canvas.paste(A2,(margin+(colw-A2.width)//2,y1)); canvas.paste(B2,(margin+colw+gap+(colw-B2.width)//2,y1)); canvas.paste(C2,((W-C2.width)//2,y2))
p=OUT/'Supplementary_Figure_S1_preview_v4.png'; canvas.save(p,dpi=(300,300)); print(p,canvas.size)
