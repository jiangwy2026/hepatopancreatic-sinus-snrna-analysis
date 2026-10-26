# Manuscript revision notes

## Result structure

The revised manuscript retains six main result sections and six main figures.

| Result section | Figure | Main content |
|---|---|---|
| 3.1 | Figure 1 | Hs is the major proliferation center after hemospasia |
| 3.2 | Figure 2 | Hs hemolymph has a distinct molecular and cellular identity |
| 3.3 | Figure 3 | Hs-plasma metabolites and betaine validation |
| 3.4 | Figure 4 | Dynamic Hs-plasma protein response |
| 3.5 | Figure 5 | Hs single-nucleus atlas, composition and PH-to-GH ordering |
| 3.6 | Figure 6 | PH substate remodeling, external progenitor concordance and cycle-independent identity |

## Integration of the former section 3.7

The standalone 3.7 heading was removed.

- The candidate niche-communication analysis was moved into section 3.5 after the composition result.
- The external LyO progenitor program and cell-cycle-adjusted identity analysis were moved into section 3.6.
- No seventh result chapter remains.

## Figure 5 revision

The source image contains panels A-G, but the previous legend stopped at F. The revised legend now covers:

- A: experimental workflow.
- B: QC metrics for 9,658 nuclei.
- C: all-nucleus t-SNE and 10 states.
- D: state marker heatmap.
- E: control and hemospasia nuclei on identical coordinates.
- F: composition and log2 fold change.
- G: PH-to-GH pseudotime ordering and marker dynamics.

Panel G is explicitly described as relative transcriptional order rather than lineage tracing.

## Figure 6 revision

The revised Figure 6 legend now matches all six panels:

- A: G1/S/G2M composition across states.
- B: PH substate composition.
- C: PH pseudotime marker heatmap.
- D: PH-only t-SNE, including the substate-centroid trajectory.
- E: external LyO progenitor and hemolymph prohemocyte UCell programs.
- F: cell-cycle-adjusted progenitor identity across Hs PH and hemolymph PH.

## Supplementary-table numbering

Existing Tables S1-S6 were retained. New tables were assigned S7-S18; Table S19 is deferred:

| Table | Content |
|---|---|
| S7 | Sample metadata |
| S8 | Nucleus-level metadata |
| S9 | QC and sequencing summary |
| S10 | Hs state markers |
| S11 | PH substate markers |
| S12 | Cell-state composition statistics |
| S13 | Control-versus-hemospasia differential expression |
| S14 | Trajectory and pseudotime source data |
| S15 | External progenitor program scoring |
| S16 | Ligand-receptor communication |
| S17 | Ligand-receptor orthology mapping |
| S18 | ITGB1 expression summary |


## Data limits recorded in the tables

- Each snRNA-seq condition is one pooled library (n=1).
- Nucleus-level P values are exploratory and are not biological-replicate estimates.
- Doublet scores are DoubletFinder pANN values; all retained nuclei were classified as Singlets.
- The `leiden_cluster` field contains the delivered vendor cluster ID because the manuscript retains vendor-derived state labels.
- Candidate communication does not establish functional signaling.
