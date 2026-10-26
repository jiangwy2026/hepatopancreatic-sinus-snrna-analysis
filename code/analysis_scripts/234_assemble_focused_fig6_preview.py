from PIL import Image,ImageOps,ImageDraw,ImageFont
import os
R=r'D:\桌面\single cell'; P=os.path.join(R,'analysis_output','figures_pub','fig5_parts'); S=os.path.join(R,'analysis_output','figures_pub','sc_summary_new_names','individual_previews'); O=os.path.join(R,'文章配图','最终_2026-10-20','主图排版预览'); os.makedirs(O,exist_ok=True)
W,H,HD=1500,1000,105
TF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',40); PF=ImageFont.truetype(r'C:\Windows\Fonts\timesbd.ttf',48); CF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',28)
panels=[('a','Fig5_qc_standardized_composition_preview.png','Composition change after QC standardization'),('b','Fig5_PH_DPT_pseudotime_density_preview.png','PH DPT pseudotime distribution'),('c','Fig5_PH_DPT_pseudotime_heatmap_preview.png','PH pseudotime-ordered marker dynamics'),('d','Fig5_R_PMlike_substate_focus_preview.png','R and PM-like substate changes'),('e','Fig_LyO_mapping_bubble_preview.png','Lymphoid-organ signature mapping'),('f',os.path.join(S,'Fig6_alluvial_supported_only_rightp_close_preview.png'),'Supported LyO programme-organ flows')]
C=Image.new('RGB',(W*2,HD+H*3),(255,255,255)); D=ImageDraw.Draw(C); title='Figure 6 | Treatment-associated state remodeling and lymphoid-organ context'; D.text(((C.width-D.textlength(title,font=TF))/2,20),title,fill=(0,0,0),font=TF)
for k,(lab,fn,cap) in enumerate(panels):
 r,c=divmod(k,2); x=c*W; y=HD+r*H; D.text((x+18,y+12),lab,fill=(0,0,0),font=PF); D.text((x+72,y+22),cap,fill=(55,55,55),font=CF)
 im=Image.open(os.path.join(P,fn) if not os.path.isabs(fn) else fn).convert('RGB'); im=ImageOps.contain(im,(W-45,H-HD-45),Image.Resampling.LANCZOS); C.paste(im,(x+(W-im.width)//2,y+HD+(H-HD-im.height)//2))
p=os.path.join(O,'Fig6_focused_main_preview.png'); C.save(p,format='PNG',dpi=(300,300)); print(p)
