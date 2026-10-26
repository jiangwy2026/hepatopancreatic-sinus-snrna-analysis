# Analysis script index

This directory is an archive of project analysis scripts. Script numbers indicate the approximate execution order; several scripts are exploratory or superseded and are retained for auditability.

## Primary entry points

| Stage | Scripts |
|---|---|
| Initial matrix loading and QC | `01_load_qc.py`, `02_markers_recluster.py`, `03_knn_merge.py`, `04_panel_and_plots.py`, `05_detection_rates.py` |
| SoupX correction and rebuild | `20a_soup_load.py` through `26_trajectory_soupx.py` |
| State markers and annotation | `274_fig2a_all_states.py`, `309_export_fig1_final_panels.py` |
| PH substates and trajectory | `248_step2_trajectory_projection.py`, `249_step2_group_stats.py`, `288_ph_substate_composition.py`, `289_ph_substate_marker_check.py`, `307_export_fig2d_ph_tsne_final.py` |
| External progenitor programs | `250_step3_lyo_program_mirror.py`, `293_cross_organ_adjusted_progenitor.py`, `295_supp14_lyo_marker_mirror.py` |
| Candidate communication | `231_step1_communication.py`, `232_step1_communication_validated.py`, `247_step1_communication_sonicparanoid.py`, `281_lr_sensitivity.py`, `286_ligand_specificity.py`, `268_supp7_receptor_group.py` |
| S7-S19 table generation | `324a_build_S7_S12.py`, `324b_build_S13_S18.py`, `324c_build_S19.py`, `325_combine_validate_tables.py`, `326_validate_S7_S19_workbook.py` |
| Manuscript rebuild and validation | `327_revise_manuscript_final.py`, `328_final_validation.py` |

## Path assumptions

Scripts retain the original local Windows path prefix used during analysis. Before running on another machine, replace the project root and input-output paths with local paths. Scripts are provided as the exact analysis record, not as a polished installable package.
