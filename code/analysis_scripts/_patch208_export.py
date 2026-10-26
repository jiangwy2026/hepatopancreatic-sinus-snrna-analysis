p=r"D:\桌面\single cell\analysis_output\scripts\208_three_column_alluvial.py"
s=open(p,encoding="utf-8").read()
if 'from PIL import Image' not in s:
 s=s.replace('from matplotlib.patches import PathPatch, Patch','from matplotlib.patches import PathPatch, Patch\nfrom PIL import Image')
old=''' p=os.path.join(PREV,"Fig6_alluvial_threecolumn_compact_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)'''
new=''' p=os.path.join(PREV,"Fig6_alluvial_threecolumn_compact_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); print("saved",p)
 deliv=os.path.join(r"D:\\桌面\\single cell\\文章配图", "最终_2026-10-19", "四器官模型_河流图"); os.makedirs(deliv,exist_ok=True)
 stem=os.path.join(deliv,"Fig6_alluvial_threecolumn_compact")
 fig.savefig(stem+".svg",bbox_inches="tight",pad_inches=0.04,facecolor="white")
 fig.savefig(stem+".tiff",dpi=1200,bbox_inches="tight",pad_inches=0.04,facecolor="white",pil_kwargs={"compression":"tiff_lzw"})
 with Image.open(stem+".tiff") as im:
  if "A" in im.getbands():
   rgba=im.convert("RGBA"); bg=Image.new("RGB",rgba.size,(255,255,255)); bg.paste(rgba,mask=rgba.getchannel("A")); tmp=stem+".tmp.tiff"; bg.save(tmp,format="TIFF",dpi=(1200,1200),compression="tiff_lzw"); os.replace(tmp,stem+".tiff")
 print("saved",stem+".svg"); print("saved",stem+".tiff")
 plt.close(fig)'''
if old not in s: raise SystemExit('old block not found')
s=s.replace(old,new)
open(p,"w",encoding="utf-8").write(s)
print("patched")
