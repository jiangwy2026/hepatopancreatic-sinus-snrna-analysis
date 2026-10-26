p=r"D:\桌面\single cell\analysis_output\scripts\190_violin_all_volcano_genes.py"
s=open(p,encoding="utf-8").read()
s=s.replace('fig = plt.figure(figsize=(7.2, 4.7), facecolor="white")','fig = plt.figure(figsize=(7.2, 4.9), facecolor="white")')
s=s.replace('left=0.10, right=0.985, top=0.86, bottom=0.11,','left=0.10, right=0.985, top=0.86, bottom=0.15,')
s=s.replace('ax.set_xticks([1, 2]); ax.set_xticklabels(["Ctrl", "Hemo"], fontsize=5.0)','ax.set_xticks([1, 2]); ax.set_xticklabels(["Control", "Hemospasia"], fontsize=4.7, rotation=20, ha="right")')
open(p,"w",encoding="utf-8").write(s)
print("patched")
