import sys, io, re
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
syms = set(str(s) for s in A.var["symbol"])
s_genes  = ["PCNA","RPA1","RPA2","RPA3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC6","CDT1","ORC1","ORC2","ORC3","ORC4","ORC5","POLA1","POLA2","POLD1","PRIM1","PRIM2","FEN1","LIG1","TYMS","RRM1","RRM2","DHFR","GINS1","GINS2","H2A","H2B","H3","H4","H1","HIRA","ASF1","CHAF1A"]
g2m_genes= ["CDK1","CDC2","CCNB","CCNB1","CCNB2","CCNB3","AURKA","AURKB","PLK1","BUB1","BUB3","MAD2","TOP2A","TPX2","KIF11","KIF15","KIF20A","KIF22","KIF23","CENPA","CENPB","SMC2","SMC4","ANLN","RACGAP1","ECT2","CDC20","CDC25","CKS1","CKS2","TUBB","TUBA","ESPL1","NDC80","NUF2","INCENP","PRC1","ASPM","MKI67","UBE2C","UBE2S","CENPE","CENPF","Pka-C1","Patronin"]
def hit(gl): return [g for g in gl if g in syms]
hs, hg = hit(s_genes), hit(g2m_genes)
print("S 期基因命中:", hs)
print("G2M 期基因命中:", hg)
open(out + r"\cellcycle_gene_hits.txt","w",encoding="utf-8").write("S: %s\nG2M: %s\n" % (hs, hg))
