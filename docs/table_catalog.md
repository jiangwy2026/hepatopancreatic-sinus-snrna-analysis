# Supplementary table catalog

Existing Tables S1-S6 are retained. The current release includes Tables S7-S18. Table S19 is deferred.

| Sheet | Rows | Columns | Title |
|---|---:|---:|---|
| S7 | 2 | 15 | Table S7. Sample metadata |
| S8 | 9,658 | 19 | Table S8. Nucleus-level metadata |
| S9 | 2 | 12 | Table S9. QC and sequencing summary |
| S10 | 2,131 | 9 | Table S10. Marker genes for Hs states |
| S11 | 226 | 8 | Table S11. PH substate marker genes |
| S12 | 10 | 13 | Table S12. Cell-state composition statistics source data |
| S13 | 1,771 | 10 | Table S13. Differential expression: control vs hemospasia |
| S14 | 40 | 10 | Table S14. Trajectory and pseudotime source data |
| S15 | 54 | 11 | Table S15. External progenitor program scoring |
| S16 | 104 | 13 | Table S16. Ligand-receptor communication |
| S17 | 21 | 6 | Table S17. LR orthology mapping |
| S18 | 14 | 6 | Table S18. ITGB1 expression summary |

## Notes

- Table S1-S6 are existing tables and were not renumbered.
- New supplementary tables start at S7.
- Each snRNA-seq condition is represented by one pooled library (n=1); nucleus-level P values are exploratory.
- DoubletFinder pANN scores were joined by Cell Ranger barcode; all retained nuclei were classified as Singlets.
- The requested leiden_cluster column contains the delivered vendor cluster ID because vendor-derived state labels were retained.

## Definition and limitation notes

### S7
- Library-level replication is n=1 per condition; no between-library variance can be estimated.

### S8
- UCell_PH1 and UCell_PH2 are eLife-derived cell-cycle and progenitor-identity program scores.
- All retained nuclei were classified as Singlets by the delivered DoubletFinder output.

### S12
- EXPLORATORY ONLY: each condition contains n=1 pooled library. Nucleus-level P values and FDR values do not estimate biological replicate variance and must not be interpreted as confirmatory group-level replication.

### S13
- Each condition is one pooled library (n=1). These nucleus-level tests are exploratory and do not substitute for biological replication.

### S15
- Scores and overlaps are descriptive and do not establish cell identity or lineage equivalence.

### S16
- subsampling_top3_frequency is a source-to-target stability metric repeated across ligand-receptor pairs sharing the same source and target.
- Candidate communication does not establish functional signaling.

### S17
- BLASTP best hits are evidence of sequence-level orthology support, not definitive one-to-one orthology.
