import pandas as pd, gzip, sys, os
sys.stdout.reconfigure(encoding="utf-8")
w = os.path.join(os.environ["TEMP"], "blastwork")
b = pd.read_csv(os.path.join(w, "new_to_old.tsv"), sep="\t", header=None,
                names=["q","s","pident","length","qlen","slen","evalue","bitscore","qcovs"])
best = b.sort_values("bitscore", ascending=False).groupby("q", as_index=False).first()
rows = []
with gzip.open(r"D:\桌面\single cell\analysis_output\GCF_003789085.1_feature_table.txt.gz","rt",encoding="utf-8",errors="ignore") as f:
    hdr = f.readline().lstrip("#").strip().split("\t")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        d = dict(zip(hdr, p))
        if d.get("feature") != "mRNA": continue
        rows.append((d.get("related_accession","").split(".")[0], d.get("GeneID",""), d.get("name","")))
ann = pd.DataFrame(rows, columns=["xp","gene_id","product"]).drop_duplicates("xp").set_index("xp")
best["old_xp"] = best.s.str.split(".").str[0]
best["gene_id"] = best.old_xp.map(lambda x: ann.loc[x,"gene_id"] if x in ann.index else "")
best["product"] = best.old_xp.map(lambda x: ann.loc[x,"product"] if x in ann.index else "")
out = best[["q","old_xp","gene_id","product","pident","qcovs","evalue"]]
out.to_csv(r"D:\桌面\single cell\analysis_output\step1\new_xp_to_geneid.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 60)
print(out.to_string(index=False))
