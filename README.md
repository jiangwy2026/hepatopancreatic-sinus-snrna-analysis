# Hepatopancreatic sinus snRNA-seq analysis

This repository contains the reproducible analysis code, supplementary source tables, and figure panel exports for a study of the hepatopancreatic sinus (Hs) as a hematopoietic niche in *Litopenaeus vannamei*.

The analysis combines BrdU/FACS validation, plasma metabolomics and proteomics, and snRNA-seq of the Hs from control and 6-h hemospasia groups. The snRNA-seq analysis resolves 10 Hs states, four PH substates, a PH-to-GH transcriptional ordering, external progenitor-program concordance, and candidate niche ligand-receptor communication.

## Repository contents

- `code/analysis_scripts/` - analysis scripts used across preprocessing, quality control, clustering, marker detection, trajectory analysis, communication, cross-organ scoring, table construction and validation.
- `tables/Supplementary_Tables_S7-S18_2026-10-26.xlsx` - combined supplementary workbook. Tables S1-S6 are existing tables and are not duplicated here.
- `tables/source_data/` - CSV versions of the source data used to build Tables S7-S18.
- `figures/Figure5_panels/` - editable panel exports for manuscript Figure 5.
- `figures/Figure6_panels/` - editable panel exports for manuscript Figure 6.
- `docs/analysis_workflow.md` - command-level analysis order and software requirements.
- `docs/table_catalog.md` - table-by-table row counts, column counts and content map.
- `docs/revision_notes.md` - manuscript and figure-legend revision record.
- `docs/data_availability.md` - public repository accessions.

## Principal analysis outputs

- 19,795 Cell Ranger nuclear barcodes; 9,658 nuclei retained after quality and doublet filtering (5,349 control; 4,309 hemospasia).
- Ten states: PH, GH, TS, MH2, F, R, B, M, MyoEpi and PM-like.
- Four PH substates: PH_s0 (ECM/ribosome-high), PH_s1 (progenitor-like), PH_s2 (immune/ECM-primed) and PH_s3 (TGM1-low).
- Figure 5 contains panels A-G, including PH-to-GH pseudotime ordering.
- Figure 6 contains panels A-F, including PH cell-cycle composition, substate composition, PH-only t-SNE, external UCell programs and cross-organ adjusted progenitor identity.

## Environment

The primary Python environment is described in `environment/requirements.txt`. The analysis used Python 3.12, Scanpy 1.12.4, AnnData 0.13.4, SoupX 1.6.2 in R, Cell Ranger 8.0.0, Seurat 4.1 and DoubletFinder 2.0.4.

## Reproducibility notes

- Raw sequencing reads, BAM files, count matrices and large AnnData objects are not redistributed in this repository. See `docs/data_availability.md`.
- The two snRNA-seq conditions are each represented by one pooled library. Nucleus-level tests are exploratory and do not estimate biological replicate variance.
- The `leiden_cluster` field in Table S8 contains the delivered vendor cluster identifier because the manuscript retains vendor-derived state labels rather than a de novo Leiden partition.
- Candidate ligand-receptor results represent expression- and orthology-based communication inference, not experimentally validated signaling.

## Data availability

- snRNA-seq BioProject: PRJNA1545411
- MetaboLights: MTBLS1478
- ProteomeXchange/iProX: PXD079762

## Citation

Please cite the accompanying manuscript when available. Repository author and license information should be completed by the corresponding author before public release.

## Code availability

All analysis scripts and supplementary source data are publicly available at https://github.com/jiangwy2026/hepatopancreatic-sinus-snrna-analysis (release `v2026.10.26-final`).
