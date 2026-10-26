from PIL import Image,ImageOps,ImageDraw,ImageFont
import os
R=r'D:\桌面\single cell'; P=os.path.join(R,'analysis_output','figures_pub','fig5_parts'); F=os.path.join(R,'analysis_output','figures_pub','final_previews'); S=os.path.join(R,'analysis_output','figures_pub','sc_summary_new_names','individual_previews'); O=os.path.join(R,'文章配图','最终_2026-10-20','主图排版预览'); os.makedirs(O,exist_ok=True)
W,H,HD=1200,850,100
TF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',38); PF=ImageFont.truetype(r'C:\Windows\Fonts\timesbd.ttf',42); CF=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',26)
F5=[('a','Fig5a_study_design_preview.png','Study design'),('b','Fig5b_qc_preview.png','QC metrics'),('c',os.path.join(F,'Fig1b_TSNE_all_groups_preview.png'),'Cell-state atlas'),('d',os.path.join(F,'Fig1c_TSNE_control_and_hemospasia_preview.png'),'Group t-SNE'),('e','Fig5_marker_bubble_10groups_preview.png','Marker identity'),('f','Fig5_dynamic_composition_cellcycle_preview.png','Composition and cell cycle'),('g','Fig5_qc_standardized_composition_preview.png','QC-standardized composition'),('h','Fig5_R_PMlike_substate_focus_preview.png','R/PM-like substates'),('i','Fig5_PH_DPT_pseudotime_density_preview.png','PH DPT pseudotime')]
F6=[('a','Fig_LyO_mapping_bubble_preview.png','LyO mapping'),('b',os.path.join(S,'Fig6_alluvial_supported_only_rightp_close_preview.png'),'Organ programme flows'),('c','Fig5_parent_program_effect_heatmap_preview.png','Program changes'),('d','Fig5_coexpression_module_heatmap_preview.png','Co-expression modules'),('e','Fig5_focused_LR_group_preview.png','Focused ligand-receptor'),('f','Fig5_DEG_sharing_upset_preview.png','Shared/unique DEGs')]
def build(panels,title,name):
 C=Image.new('RGB',(W*3,HD+H*3),(255,255,255)); D=ImageDraw.Draw(C); D.text(((C.width-D.textlength(title,font=TF))/2,20),title,fill=(0,0,0),font=TF)
 for k,(lab,fn,cap) in enumerate(panels):
  r,c=divmod(k,3); x=c*W; y=HD+r*H; D.text((x+16,y+12),lab,fill=(0,0,0),font=PF); D.text((x+64,y+22),cap,fill=(55,55,55),font=CF)
  if fn is None: continue
  im=Image.open(os.path.join(P,fn) if not os.path.isabs(fn) else fn).convert('RGB')
  if os.path.basename(fn) in ('Fig5a_study_design_preview.png','Fig5b_qc_preview.png'):
   ImageDraw.Draw(im).rectangle((0,0,85,85),fill='white')
  im=ImageOps.contain(im,(W-40,H-HD-40),Image.Resampling.LANCZOS); C.paste(im,(x+(W-im.width)//2,y+HD+(H-HD-im.height)//2))
 p=os.path.join(O,name); C.save(p,format='PNG',dpi=(300,300)); return C
a=build(F5,'Figure 5 | Single-cell landscape and treatment-associated remodeling','Fig5_main_preview.png')
b=build(F6,'Figure 6 | Lymphoid-organ context and candidate hematopoietic axis','Fig6_main_preview.png')
G=Image.new('RGB',(a.width,a.height+b.height+50),(255,255,255)); G.paste(a,(0,0)); G.paste(b,(0,a.height+50)); G.save(os.path.join(O,'Main_Figures_overview_preview.png'),format='PNG',dpi=(300,300)); print('saved')
