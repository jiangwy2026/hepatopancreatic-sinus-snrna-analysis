# -*- coding: utf-8 -*-
"""267_supp6_prohemocyte.py — PH vs eLife prohemocyte (PH1/PH2) direct comparison."""
import os, collections
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Supplementary_niche")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})
def save(fig,n):
    fig.savefig(os.path.join(DEST,"svg",n+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",n+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",n+".png"),dpi=220,bbox_inches="tight"); plt.close(fig); print("saved",n,flush=True)
def letter(ax,s): ax.text(-0.10,1.05,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

a=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
vn=set(map(str,a.var_names)); sym=collections.defaultdict(list)
for i,s in enumerate(a.var["sym"].astype(str)): sym[s.lower()].append(a.var_names[i])
def map_set(df,n=50):
    d=df.sort_values("Log2FC",ascending=False).head(n)
    genes=[]; labels=[]
    for x in d["Gene Name"].astype(str):
        v=x if x in vn else (sym[x.lower()][0] if len(sym.get(x.lower(),[]))==1 else None)
        if v and v not in genes: genes.append(v); labels.append(x)
    return genes,labels
p1,_=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp3-v2.csv")),50)
p2,_=map_set(pd.read_csv(os.path.join(OUT,"elife_ref","elife-80127-supp6-v2.csv")),50)
print("eLife PH1 mapped",len(p1),"PH2 mapped",len(p2))

def ucell(genes,chunk=400):
    idx=[list(a.var_names).index(g) for g in genes]; k=len(idx); N=a.n_vars; out=np.zeros(a.n_obs,dtype=np.float32)
    for st in range(0,a.n_obs,chunk):
        en=min(a.n_obs,st+chunk); X=np.asarray(a[st:en].X.todense(),dtype=np.float32)
        rank=np.argsort(-X,axis=1,kind="stable"); r=np.empty_like(rank); r[np.arange(en-st)[:,None],rank]=np.arange(1,N+1)
        S=r[:,idx].sum(1); out[st:en]=1-(S-k*(k+1)/2)/(k*(N-k))
    return out
a.obs["ucell_eLifePH1"]=ucell(p1); a.obs["ucell_eLifePH2"]=ucell(p2)

# hypergeometric
mk=pd.read_csv(os.path.join(OUT,"markers_10groups.csv"))
def top(grp,n=50):
    g=mk[mk.grp==grp].sort_values("scores",ascending=False).head(n)
    return set(x for x in g["names"].astype(str) if x in vn)
ph=top("PH"); niche=set().union(*[top(g) for g in ["TS","MH2","MyoEpi","M"]])
rows=[]
for name,S in [("eLife PH1",set(p1)),("eLife PH2",set(p2))]:
    for lab,B in [("PH markers",ph),("niche markers",niche)]:
        k=len(S&B); M=len(vn); n=len(S); N=len(B); p=hypergeom.sf(k-1,M,n,N) if k else 1.0
        rows.append({"set":name,"compared_to":lab,"overlap":k,"expected":n*N/M,"fold":k/(n*N/M),"p":p,"genes":";".join(sorted(S&B))})
ov=pd.DataFrame(rows); ov.to_csv(os.path.join(DEST,"Supp6_hypergeom.csv"),index=False,encoding="utf-8-sig")

fig=plt.figure(figsize=(7.2,4.6)); gs=fig.add_gridspec(1,3,width_ratios=[1.5,1.0,1.0],wspace=0.45,left=0.08,right=0.97,top=0.9,bottom=0.22)
# a boxplots
ax=fig.add_subplot(gs[0])
data=[]; labels=[]; cols=[]
for key,lab,col in [("ucell_eLifePH1","eLife PH1","#4C78A8"),("ucell_eLifePH2","eLife PH2","#E45756"),("ucell_LyO_Progenitor","LyO Progenitor","#54A24B")]:
    for ct in order:
        data.append(a.obs.loc[a.obs.cell_type==ct,key].values); labels.append(ct); cols.append(col)
bp=ax.boxplot(data,patch_artist=True,showfliers=False,widths=0.62,medianprops=dict(color="#222",lw=0.6))
for p,c in zip(bp["boxes"],cols): p.set_facecolor(c); p.set_alpha(0.55); p.set_edgecolor("#444"); p.set_linewidth(0.4)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels,rotation=90,fontsize=5.2); ax.set_ylabel("UCell score")
ax.set_title("Prohemocyte programmes"); 
for i in range(1,len(order)): ax.axvline(i*3-0.5,color="#DDD",lw=0.5)
ax.legend(handles=[Patch(color="#4C78A8",label="eLife PH1"),Patch(color="#E45756",label="eLife PH2"),Patch(color="#54A24B",label="LyO Progenitor")],frameon=False,fontsize=5.2,loc="upper right")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"a")
# b hypergeom
ax=fig.add_subplot(gs[1])
sub=ov.copy(); sub["lab"]=sub.set+"\nvs "+sub.compared_to
y=np.arange(len(sub)); ax.barh(y,sub.fold,color=["#C44E52" if "PH markers" in c else "#999" for c in sub.compared_to])
ax.set_yticks(y); ax.set_yticklabels(sub.lab,fontsize=5.2); ax.invert_yaxis()
for i,r in sub.iterrows(): ax.text(r.fold+0.2,i,f"{int(r.overlap)} genes\nP={r.p:.1e}",va="center",fontsize=4.8)
ax.set_xlabel("Fold enrichment"); ax.set_title("Overlap with Hs markers")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"b")
# c PH substates along PH->GH projected pseudotime
ax=fig.add_subplot(gs[2])
t=ad.read_h5ad(os.path.join(OUT,"step2","step2_integrated.h5ad"),backed="r")
phm=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv")); phm.columns=["cell","group","subcluster","u1","u2"]
sub=phm.set_index("cell")["subcluster"]
d=t.obs.loc[t.obs.cell_type=="PH",["projected_pseudotime"]].copy(); d["sub"]=sub.reindex(d.index)
order_s=[1,0,2,3]
for pos,s_ in enumerate(order_s):
    v=d.loc[d["sub"]==s_,"projected_pseudotime"].dropna().values
    ax.boxplot([v],positions=[pos],widths=0.5,patch_artist=True,
               boxprops=dict(facecolor="#B2182B",alpha=0.6),showfliers=False,medianprops=dict(color="#222",lw=0.6))
ax.set_xticks(range(4)); ax.set_xticklabels([f"s{s_}" for s_ in order_s]); ax.set_ylabel("Projected pseudotime (PH->GH)")
ax.set_title("PH substates along trajectory")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"c")
fig.savefig(os.path.join(DEST,"png","Supp6_prohemocyte_comparison.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(DEST,"svg","Supp6_prohemocyte_comparison.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"tiff","Supp6_prohemocyte_comparison.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
print("done")
