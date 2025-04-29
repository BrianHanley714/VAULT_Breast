---
title: "Breast Cancer Analysis Report"
author: "Greg"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 2
    number_sections: true
    theme: flatly
    highlight: zenburn
---

<p align="center">
  <img src="assets/Logo_VAULT.png" alt="VAULT Logo" width="200"/>
</p>

# VAULT_Breast
## Content of this repository
This is a repository containing the data and code required to reproduce the results and figures of the manuscript "*Capturing actionable breast cancer evolution in discarded tissue*".
Below an indect how how the scripts relate to the different figures. 

## Code for reproduction of the manuscript figures
The table below links which code reproduces each manuscript figure

| Figure type   | Figure number  | Code location |
| ------------- | -------------- | --------------| 
| Main          | 1a              | None, this was generated using graphical design software|
| Main          | 1b              | [analysis/scripts/circos_cohort_overview.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 1c              | None, this was generated using graphical design softwareand QuPath|
| Main          | 1d              | None, this was generated using graphical design software|
| Main          | 1e              | [analysis/scripts/prop_leftover_tissue.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2a              | [analysis/scripts/density_histogram_rw_vault.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2b              | [analysis/scripts/assess_limitofdetection.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2c              | [analysis/scripts/density_histogram_simulated.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2d              | [analysis/scripts/density_histogram_rw_vault.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2e              | [analysis/scripts/variant_clonality_variant_count.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3a              | [analysis/scripts/Prop_cohort_w_drivers.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3b              | [analysis/scripts/biomarker_status_by_clonality.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3c              | [analysis/scripts/vault_biomarkers_clonality_genes.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3d              | [analysis/scripts/PIK3CA_lollipop.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3e              | [analysis/scripts/variant_clonality_actionability.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3f              | [analysis/scripts/driver_CCF_by_study.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3g              | [analysis/scripts/proportion_patients_w_driver_by_study.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3h              | [analysis/scripts/known_v_predicted_drivers.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3i              | [analysis/scripts/MYEOV_PTH2_lollipops.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 3j              | [analysis/scripts/clonality_selected_drivers.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|




