p=r"D:\桌面\single cell\analysis_output\scripts\190_violin_all_volcano_genes.py"
s=open(p,encoding="utf-8").read()
s=s.replace('''        det = [100 * (data[0] > 0).mean(), 100 * (data[1] > 0).mean()]
        ax.text(1.5, -0.035, "detected %.0f%% → %.0f%%" % (det[0], det[1]),
                ha="center", va="top", fontsize=4.7, color="#555555")
''','')
open(p,"w",encoding="utf-8").write(s)
print("patched")
