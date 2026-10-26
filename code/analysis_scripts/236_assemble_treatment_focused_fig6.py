from PIL import Image,ImageOps,ImageDraw,ImageFont
import os
R=r'D:\桌面\single cell'; P=os.path.join(R,'analysis_output','figures_pub','fig5_parts'); S=os.path.join(R,'analysis_output','figures_pub','sc_summary_new_names','individual_previews'); O=os.path.join(R,'文章配图','最终_2026-10-20','主图排版预览'); os.makedirs(O,exist_ok=True)
W,H,HD=1500,1000,105
TF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',40); PF=ImageFont.truetype(r'C:\Windows\Fonts\timesbd.ttf',48); CF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',28)
top=[('a','Fig5_PH_DPT_pseudotime_density_preview.png','PH DPT pseudotime distribution'),('b','Fig5_PH_DPT_pseudotime_heatmap_preview.png','PH pseudotime-ordered markers'),('c','Fig5_cellcycle_dynamic_preview.png','Cell-cycle phase by group'),('d','Fig5_parent_program_effect_heatmap_preview.png','QC-adjusted program changes'),('e','Fig5_DEG_sharing_upset_preview.png','Shared and unique DEGs'),('f','Fig5_focused_LR_group_preview.png','Focused ligand-receptor changes')]
C=Image.new('RGB',(W*3,HD+H*2+H+60),(255,255,255)); D=ImageDraw.Draw(C); title='Figure 6 | Hemospasia-induced hemocyte state remodeling'; D.text(((C.width-D.textlength(title,font=TF))/2,20),title,fill=(0,0,0),font=TF)
for k,(lab,fn,cap) in enumerate(top):
 r,c=divmod(k,3); x=c*W; y=HD+r*H; D.text((x+18,y+12),lab,fill=(0,0,0),font=PF); D.text((x+72,y+22),cap,fill=(55,55,55),font=CF)
 im=Image.open(os.path.join(P,fn)).convert('RGB'); im=ImageOps.contain(im,(W-45,H-HD-45),Image.Resampling.LANCZOS); C.paste(im,(x+(W-im.width)//2,y+HD+(H-HD-im.height)//2))
y=HD+2*H; D.text((18,y+12),'g',fill=(0,0,0),font=PF); D.text((72,y+22),'Supported LyO programme-organ flows: bridge to the four-organ model',fill=(55,55,55),font=CF)
im=Image.open(os.path.join(S,'Fig6_alluvial_supported_only_rightp_close_preview.png')).convert('RGB'); box=(C.width-50,H+10); im=ImageOps.contain(im,box,Image.Resampling.LANCZOS); C.paste(im,((C.width-im.width)//2,y+HD+(-10)))
p=os.path.join(O,'Fig6_treatment_focused_preview.png'); C.save(p,format='PNG',dpi=(300,300)); print(p)
