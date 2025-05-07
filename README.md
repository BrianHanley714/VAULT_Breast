|                            |                                  |
|----------------------------|----------------------------------|
| ![](assets/Logo_VAULT.png) | ![](assets/VAULT_Breastlogo.png) |

# VAULT_Breast

## Overview

This is a repository containing the data and code required to reproduce the results and figures of the manuscript "*Capturing actionable breast cancer evolution in discarded tissue*". Below is a summary table of how the scripts relate to the different figures.

------------------------------------------------------------------------

## Data

The datasets used in the study are available in the [data/](data/) directory.\
Please refer to the data license terms before use.

------------------------------------------------------------------------

## Code for reproduction of the manuscript figures

The table below links which code reproduces each manuscript figure

| Figure location | Figure number | Code location                                                                                                                                                                 |
|-----------------------|-------------------------|-------------------------|
| Main            | 1a            | None, this was generated using graphical design software                                                                                                                      |
| Main            | 1b            | [analysis/scripts/circos_cohort_overview.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/circos_cohort_overview.R)                               |
| Main            | 1c            | None, this was generated using graphical design software and QuPath                                                                                                           |
| Main            | 1d            | None, this was generated using graphical design software                                                                                                                      |
| Main            | 1e            | [analysis/scripts/prop_leftover_tissue.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/prop_leftover_tissue.R)                                   |
| Main            | 2a            | [analysis/scripts/density_histogram_rw_vault.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/density_histogram_rw_vault.R)                       |
| Main            | 2b            | [analysis/scripts/assess_limitofdetection.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/assess_limitofdetection.R)                             |
| Main            | 2c            | [analysis/scripts/density_histogram_simulated.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/density_histogram_simulated.R)                     |
| Main            | 2d            | [analysis/scripts/density_histogram_rw_vault.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/density_histogram_rw_vault.R)                       |
| Main            | 2e            | [analysis/scripts/variant_clonality_variant_count.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/variant_clonality_variant_count.R)             |
| Main            | 3a            | [analysis/scripts/Prop_cohort_w_drivers.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/Prop_cohort_w_drivers.R)                                 |
| Main            | 3b            | [analysis/scripts/biomarker_status_by_clonality.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/biomarker_status_by_clonality.R)                 |
| Main            | 3c            | [analysis/scripts/vault_biomarkers_clonality_genes.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/vault_biomarkers_clonality_genes.R)           |
| Main            | 3d            | [analysis/scripts/PIK3CA_lollipop.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/PIK3CA_lollipop.R)                                             |
| Main            | 3e            | [analysis/scripts/variant_clonality_actionability.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/umap.R)                                        |
| Main            | 3f            | [analysis/scripts/driver_CCF_by_study.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/variant_clonality_actionability.R)                         |
| Main            | 3g            | [analysis/scripts/proportion_patients_w_driver_by_study.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/proportion_patients_w_driver_by_study.R) |
| Main            | 3h            | [analysis/scripts/known_v_predicted_drivers.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/known_v_predicted_drivers.R)                         |
| Main            | 3i            | [analysis/scripts/MYEOV_PTH2_lollipops.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/MYEOV_PTH2_lollipops.R)                                   |
| Main            | 3j            | [analysis/scripts/clonality_selected_drivers.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/clonality_selected_drivers.R)                       |
| Main            | 4a            | [analysis/scripts/draw_trees_for_VAULT.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/draw_trees_for_VAULT.R)                                   |
| Main            | 4b            | [analysis/scripts/draw_trees_for_VAULT.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/draw_trees_for_VAULT.R)                                   |
| Main            | 4c            | [analysis/scripts/draw_trees_for_VAULT.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/draw_trees_for_VAULT.R)                                   |
| Main            | 5a            | [analysis/scripts/ki67_heterogeneity.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/ki67_heterogeneity.R)                                       |
| Main            | 5b            | [analysis/scripts/HF007_ki67_heterogeneity.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/HF007_ki67_heterogeneity.R)                           |
| Main            | 5c            | [analysis/scripts/plot_phenophylogenies.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/plot_phenophylogenies.R)                                 |
| Main            | 5d            | [analysis/scripts/HF299_phenophylogeny.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/HF299_phenophylogeny.R)                                   |
| Extended        | 1             | [analysis/scripts/Consort_Diagram.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/Consort_Diagram.R)                                             |
| Extended        | 2a,b          | [analysis/scripts/secondary_endpoints.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/secondary_endpoints.R)                                     |
| Extended        | 3a            | None, this was generated using graphical design software)                                                                                                                     |
| Extended        | 3b-p          | [analysis/scripts/secondary_endpoints.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/secondary_endpoints.R)                                     |
| Extended        | 4             | [analysis/scripts/get_matching_data.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/get_matching_data.R)                                         |
| Extended        | 5             | [analysis/scripts/neoadj_Rx_increases_ITH.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/neoadj_Rx_increases_ITH.R)                             |
| Extended        | 6             | [analysis/scripts/compare_ffpesig.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/compare_ffpesig.R)                                             |
| Extended        | 7             | [analysis/scripts/lvi_captures_met.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/lvi_captures_met.R)                                           |
| Extended        | 8             | [analysis/scripts/draw_trees_clonal_illusion.R](https://github.com/BrianHanley714/VAULT_Breast/tree/main/analysis/scripts/draw_trees_clonal_illusion.R)                       |



------------------------------------------------------------------------

## Citation

The manuscript is in review and citation will be updated upon publication

------------------------------------------------------------------------

## Contact

For questions, feedback, or collaborations, please contact:

-   Dr Brian Hanley (Bioinformatics Lead) \| [brian.hanley\@crick.ac.uk](mailto:brian.hanley@crick.ac.uk)
-   Prof Samra Turajlic (Principal Investigator) \| [samra.turajlic\@crick.ac.uk](mailto:samra.turajlic@crick.ac.uk)
-   Project Website \| [crick.ac.uk/RepSamp](https://www.crick.ac.uk/research/labs/samra-turajlic/areas-of-interest/representative-sampling)
-   Clinical Trial : [clinicaltrials.gov/VAULT](https://clinicaltrials.gov/study/NCT03832062?term=NCT03832062&rank=1)
