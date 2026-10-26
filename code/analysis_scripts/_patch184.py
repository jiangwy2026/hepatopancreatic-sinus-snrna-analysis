p=r"D:\桌面\single cell\analysis_output\scripts\184_ph_tf_screen.py"
s=open(p,encoding="utf-8").read()
s=s.replace('pd.read_csv(os.path.join(ROOT, "all.genes.marker.stat.xls.xls"), sep="\\t", low_memory=False)', 'pd.read_csv(os.path.join(ROOT, "all.genes.marker.stat.xls.xls"), sep="\\t", low_memory=False, encoding="gb18030")')
open(p,"w",encoding="utf-8").write(s)
print("patched")
