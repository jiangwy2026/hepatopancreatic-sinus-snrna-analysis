import os
import numpy as np, pandas as pd, anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
a=ad.read_h5ad(os.path.join(OUT,"step3","step3_lyo_scores.h5ad"),backed="r")
g2v={str(g).replace("ncbi_",""): g for g in a.var_names}
def expr(ids): return np.asarray(a[:,[g2v[i] for i in ids if i in g2v]].layers["lognorm"].todense()).mean(1)
vwf=expr(["113824171"]); itgb1=expr(["113830322","113805715"])
ct=a.obs.cell_type.astype(str).values; grp=a.obs.group.astype(str).values
rng=np.random.default_rng(1)
def boot(v,n=2000):
    out=[np.mean(rng.choice(v,len(v),replace=True)) for _ in range(n)]
    return np.percentile(out,[2.5,97.5])
myo=vwf[ct=="MyoEpi"]; ph=itgb1[ct=="PH"]; gh=itgb1[ct=="GH"]
print("MyoEpi n",len(myo),"VWF mean",myo.mean(),"CI",boot(myo),"pct+",(myo>0).mean()*100)
# communication bootstrap
for tgt_name,tgt in [("PH",ph),("GH",gh)]:
    vals=[]
    for _ in range(2000):
        m=rng.choice(myo,len(myo),replace=True); t=rng.choice(tgt,len(tgt),replace=True)
        vals.append(m.mean()*t.mean())
    print(f"MyoEpi->{tgt_name} prob CI",np.percentile(vals,[2.5,50,97.5]))
# group difference MyoEpi
mc=vwf[(ct=="MyoEpi")&(grp=="blank")]; mt=vwf[(ct=="MyoEpi")&(grp=="treated")]
print("MyoEpi control n",len(mc),"treated n",len(mt),"means",mc.mean(),mt.mean())
vals=[]
for _ in range(2000):
    vals.append(rng.choice(mt,len(mt),replace=True).mean()-rng.choice(mc,len(mc),replace=True).mean())
print("VWF diff CI",np.percentile(vals,[2.5,50,97.5]),"P(diff<=0)",np.mean(np.array(vals)<=0))
