p=r"D:\桌面\single cell\analysis_output\scripts\185_scenic_prune_annotation_test.py"
s=open(p,encoding="utf-8").read()
s=s.replace('adj.columns=["TF","target","importance"]','adj.columns=["TF","target","weight"]')
open(p,"w",encoding="utf-8").write(s)
print("patched")
