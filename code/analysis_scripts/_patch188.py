p=r"D:\桌面\single cell\analysis_output\scripts\188_ccnd2_TGM_mechanism_preview.py"
s=open(p,encoding="utf-8").read()
s=s.replace('ax.add_patch(FancyBboxPatch((0.18, 0.69), 0.44, 0.12, boxstyle="round,pad=0.01,rounding_size=0.02",', 'ax.add_patch(FancyBboxPatch((0.18, 0.82), 0.44, 0.09, boxstyle="round,pad=0.01,rounding_size=0.02",')
s=s.replace('ax.text(0.40, 0.75, "Candidate upstream: br / Trl / Cfp1 / Taf7 / Kdm4B",', 'ax.text(0.40, 0.865, "Candidate upstream: br / Trl / Cfp1 / Taf7 / Kdm4B",')
open(p,"w",encoding="utf-8").write(s)
print("patched")
