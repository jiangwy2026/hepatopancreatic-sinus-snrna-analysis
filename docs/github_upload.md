# GitHub upload

## Current local status

- Repository: `D:\桌面\single cell\github_release_2026-10-26`
- Branch: `main`
- Tag: `v2026.10.26`
- Current release: Tables S7-S18; Table S19 is deferred and is not included.
- Repository size: about 30 MB.

## Method 1: GitHub CLI

Install GitHub CLI if it is not already available, then authenticate:

```powershell
gh auth login
```

Create a private repository and push the prepared history:

```powershell
cd "D:\桌面\single cell\github_release_2026-10-26"
gh repo create hepatopancreatic-sinus-snrna-analysis --private --source=. --remote=origin --push
git push origin v2026.10.26
```

Use `--public` instead of `--private` only after the corresponding author has approved public release.

## Method 2: Create an empty repository in the browser

Create a repository without a README, license or `.gitignore`, then run:

```powershell
cd "D:\桌面\single cell\github_release_2026-10-26"
git remote add origin https://github.com/<USERNAME>/hepatopancreatic-sinus-snrna-analysis.git
git push -u origin main
git push origin v2026.10.26
```

## Method 3: Restore from the Git bundle

```powershell
git clone "D:\桌面\single cell\修订交付_2026-10-26_v2\hepatopancreatic_sinus_snrna_analysis_S7_S18_final.bundle" hepatopancreatic-sinus-snrna-analysis
cd hepatopancreatic-sinus-snrna-analysis
git remote add origin https://github.com/<USERNAME>/hepatopancreatic-sinus-snrna-analysis.git
git push -u origin main
git push origin v2026.10.26
```

## What will be uploaded

- Analysis scripts
- Tables S7-S18 combined workbook
- Tables S7-S18 source CSV files
- Figure 5 and Figure 6 PNG/SVG/PDF panels
- Analysis workflow, table catalog, data availability and validation records

The manuscript DOCX, raw FASTQ/BAM files, H5AD objects and large TIFF files are intentionally excluded.
