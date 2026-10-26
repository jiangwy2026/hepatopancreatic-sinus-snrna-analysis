# GitHub upload

The local repository is initialized on branch `main`, tagged `v2026.10.26`, and contains no files larger than 5 MB.

## Current status

- Local repository: `D:\桌面\single cell\github_release_2026-10-26`
- Git bundle: `D:\桌面\single cell\修订交付_2026-10-26_v2\hepatopancreatic_sinus_snrna_analysis_2026-10-26.bundle`
- Latest commit before this document: `2468689`
- Tag: `v2026.10.26`
- Visibility: not yet published

## Option A: GitHub CLI after login

```powershell
cd "D:\桌面\single cell\github_release_2026-10-26"
gh auth login
gh repo create hepatopancreatic-sinus-snrna-analysis --private --source=. --remote=origin --push
git push origin v2026.10.26
```

## Option B: create an empty GitHub repository in the browser

1. Create an empty repository named `hepatopancreatic-sinus-snrna-analysis`.
2. Do not initialize it with a README, license or `.gitignore`.
3. Run:

```powershell
cd "D:\桌面\single cell\github_release_2026-10-26"
git remote add origin https://github.com/<USERNAME>/hepatopancreatic-sinus-snrna-analysis.git
git push -u origin main
git push origin v2026.10.26
```

## Option C: restore from the bundle

```powershell
git clone hepatopancreatic_sinus_snrna_analysis_2026-10-26.bundle hepatopancreatic-sinus-snrna-analysis
cd hepatopancreatic-sinus-snrna-analysis
git remote add origin https://github.com/<USERNAME>/hepatopancreatic-sinus-snrna-analysis.git
git push -u origin main
git push origin v2026.10.26
```

Do not upload the manuscript DOCX, raw FASTQ/BAM files, H5AD objects or large TIFF files to this code repository unless the corresponding author has confirmed the journal and institutional data-sharing requirements.
