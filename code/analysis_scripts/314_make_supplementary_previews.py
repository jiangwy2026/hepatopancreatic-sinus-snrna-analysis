# -*- coding: utf-8 -*-
"""Compose supplementary figure previews without altering source panels."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(r'D:\桌面\single cell')
SUP=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary'
OUT=ROOT/'analysis_output'/'figures_pub'/'2026-10-23_Fig1_Fig2_final'/'Supplementary_previews'
OUT.mkdir(parents=True,exist_ok=True)
font_path=r'C:\Windows\Fonts\times.ttf'; font_bold=r'C:\Windows\Fonts\timesbd.ttf'

def load(p):
    im=Image.open(p).convert('RGBA')
    bg=Image.new('RGBA',im.size,(255,255,255,255)); bg.alpha_composite(im)
    return bg.convert('RGB')

def contain(im,w,h):
    im=im.copy(); im.thumbnail((w,h),Image.Resampling.LANCZOS); return im

def make_s1():
    margin=40; gap=30; header=105; W=2400
    A=load(SUP/'Supp_PH_substate_markers.png'); B=load(SUP/'Supp_PH_identity_vs_cellcycle.png'); C=load(SUP/'Supp_LyO_marker_mirror.png')
    top_h=1000; colw=(W-2*margin-gap)//2
    A2=contain(A,colw,top_h); B2=contain(B,colw,top_h)
    C2=contain(C,W-2*margin,1500)
    row1_y=margin+header
    bottom_y=row1_y+top_h+gap
    H=bottom_y+C2.height+margin
    canvas=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(canvas)
    d.text((W//2,35),'Supplementary Figure S1 | PH substate identity evidence',font=ImageFont.truetype(font_bold,38),fill='#0f172a',anchor='mm')
    # top labels and images
    xA=margin+(colw-A2.width)//2; xB=margin+colw+gap+(colw-B2.width)//2
    d.text((margin+8,row1_y-2),'A',font=ImageFont.truetype(font_bold,32),fill='black')
    d.text((margin+colw+gap+8,row1_y-2),'B',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(A2,(xA,row1_y+38)); canvas.paste(B2,(xB,row1_y+38))
    d.text((margin+8,bottom_y-2),'C',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(C2,((W-C2.width)//2,bottom_y+38))
    p=OUT/'Supplementary_Figure_S1_preview.png'; canvas.save(p,dpi=(300,300)); print(p,canvas.size)

def make_s2():
    margin=40; gap=28; header=105; W=2400
    A=load(SUP/'Supp_Niche_marker_dotplot.png'); B=load(SUP/'Supp_Niche_LR_pairs.png'); C=load(SUP/'Supp_Niche_orthogroup_presence.png'); D=load(SUP/'Supp_Niche_receptor_expression.png'); E=load(SUP/'Supp_Niche_communication_sensitivity.png')
    colw=(W-2*margin-gap)//2; top_h=720
    A2=contain(A,colw,top_h); B2=contain(B,colw,top_h)
    mid_h=1000
    C2=contain(C,320,mid_h); D2=contain(D,W-2*margin-320-gap,mid_h)
    E2=contain(E,W-2*margin,950)
    y1=margin+header; y2=y1+top_h+gap; y3=y2+mid_h+gap
    H=y3+E2.height+margin
    canvas=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(canvas)
    d.text((W//2,35),'Supplementary Figure S2 | Candidate Hs niche signaling',font=ImageFont.truetype(font_bold,38),fill='#0f172a',anchor='mm')
    xA=margin+(colw-A2.width)//2; xB=margin+colw+gap+(colw-B2.width)//2
    d.text((margin+8,y1-2),'A',font=ImageFont.truetype(font_bold,32),fill='black')
    d.text((margin+colw+gap+8,y1-2),'B',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(A2,(xA,y1+38)); canvas.paste(B2,(xB,y1+38))
    # mid row: C narrow, D wide, centered as a group
    groupw=C2.width+gap+D2.width; x0=(W-groupw)//2
    d.text((x0-28,y2-2),'C',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(C2,(x0,y2+38))
    d.text((x0+C2.width+gap+8,y2-2),'D',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(D2,(x0+C2.width+gap,y2+38))
    d.text((margin+8,y3-2),'E',font=ImageFont.truetype(font_bold,32),fill='black')
    canvas.paste(E2,((W-E2.width)//2,y3+38))
    p=OUT/'Supplementary_Figure_S2_niche_preview.png'; canvas.save(p,dpi=(300,300)); print(p,canvas.size)

make_s1(); make_s2()
