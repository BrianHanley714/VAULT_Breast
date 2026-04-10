rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(ComplexHeatmap)
library(qs)
library(openxlsx)
library(tidyverse)
library(qs)

# LOAD DATA ---------------------------------------------------------------
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
TREES_PATHS = file.path(BASE, "data", "Trees")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")
FLOW_PATH = file.path(BASE, "data", "fixed_flow", "fixed_flow.txt")
VAULT_KI67_SORTED = file.path(BASE, "data", "metadata", "Ki67_sorted_cases.xlsx")
KI67_LOC =  file.path(BASE, "data", "image_analysis", "Ki67_pathologist_scores.txt")
ENDPOINT_DATA = file.path(BASE, "data", "metadata", "VAULT_endpoint_data.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))
source(file.path(BASE, "src", "custom_tree_plotting.R"))
source(file.path(BASE, "src", "conipher_helper_functions.R"))


# LOAD DATA ---------------------------------------------------------------


flow_data = read.delim(FLOW_PATH)
clinical_data = read.delim(VAULT_CLIN_PATH)
Ki67cases = read.xlsx(VAULT_KI67_SORTED, colNames = F)[,1]
all_variants = qread(VAULT_VARIANTS_PATH)
mat_out = read.delim(KI67_LOC)
endpoint_data= read.delim(ENDPOINT_DATA)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]

# Organise data -----------------------------------------------------------



mat_out
flow_data = flow_data%>%dplyr::mutate(Trial.ID = substr(Trial_id, 1, 5))
names(flow_data)[names(flow_data) %in% names(clinical_data)]
clinical_data$ER_alred = as.numeric(substr(clinical_data$ER_Result, 1, 1))
clinical_data$PR_alred = as.numeric(substr(clinical_data$PR_Result, 1, 1))

clinical_data = left_join(clinical_data, flow_data, by = "Trial.ID")


clinical_data = clinical_data%>%
  #dplyr::filter(Trial.ID %in% substr(ck_included_cases_vault, 1, 5))%>%
  mutate(ER_low = ER_alred<=3,
         ER_index = ER/CK8.18,
         PR_low = PR_alred<=3,
         PR_index = PR/CK8.18,
         HER2_index = HER2/CK8.18,
         Her.2_binned = if_else(Her.2_Result == "3+", "Positive", "Negative/borderline")
  )

clinical_data = clinical_data%>%
  dplyr::filter(Trial.ID %in% substr(Ki67cases, 1, 5))


sub_clonal_count = all_variants%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                clonality == "SUBCLONAL")%>%
  dplyr::count(Tumor_Sample_Barcode)%>%
  mutate(Trial_ID = substr(Tumor_Sample_Barcode, 1, 5))

total_count = all_variants%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
  dplyr::count(Tumor_Sample_Barcode)%>%
  mutate(Trial_ID = substr(Tumor_Sample_Barcode, 1, 5))

driver_count = all_variants%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                Consensus_annotated != "Non_Driver")%>%
  dplyr::count(Tumor_Sample_Barcode)%>%
  mutate(Trial_ID = substr(Tumor_Sample_Barcode, 1, 5))

order = 
  mat_out%>%
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         order = sum(order))%>%
  arrange(Trial_ID, order)%>%
  dplyr::filter(!duplicated(filenames))%>%
  ungroup()%>%pull(filenames)

mat_out = mat_out%>%  
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         region_mean = sum(order))%>%
  group_by(Trial_ID)%>%
  mutate(max = max(region_mean),
         min = min(region_mean),
         mean = mean(region_mean),
         range = max-min)
mat_out = mat_out%>%
  mutate(Trial.ID = Trial_ID)%>%
  group_by(Trial.ID)%>%
  reframe(max_ki67 = max(max),
          min_ki67 = min(min),
          mean_ki67 = mean(mean)
  )

clinical_data = left_join(clinical_data, mat_out, by = "Trial.ID")

clinical_data$total_count = total_count$n[match(clinical_data$Trial.ID, total_count$Trial_ID)]
clinical_data$driver_count = driver_count$n[match(clinical_data$Trial.ID, driver_count$Trial_ID)]
clinical_data$subclonal_mut_count = sub_clonal_count$n[match(clinical_data$Trial.ID, sub_clonal_count$Trial_ID)]
ki67_cases_driver_subclones = c("HF299", "HF039", "HF414", "HF105", "HF168")

clinical_data$Tumour_Weight = endpoint_data$Tumour.leftover.leftover.tissue.weight..grams.[match(clinical_data$Trial.ID, endpoint_data$Trial.ID)]


clinical_data = clinical_data%>%
  mutate(driver_subclone_identified = Trial.ID %in% ki67_cases_driver_subclones)

clinical_data = clinical_data%>%
  mutate(Ki67_index = KI67.ALL/CK8.18)%>%
  mutate(range_ki67 = max_ki67 - min_ki67)




comparators = c("Invasive.tumour.size..mm.","Ki67_index", "Tumour_Weight", "total_count", "Pathology.Grading.Components_M",
                "driver_count",
                "max_ki67", 
                "min_ki67", 
                "range_ki67",
                "mean_ki67",
                "subclonal_mut_count",
                "TNM8_Stage")
plotting_mat = clinical_data%>%dplyr::select(comparators)%>%as.matrix()

check_stats = function(return = "pval"){sapply(comparators, function(comparator){
  wc_test = wilcox.test(clinical_data%>%
                          dplyr::filter(driver_subclone_identified)%>%
                          pull(comparator)%>%as.numeric()
                        ,
                        clinical_data%>%
                          dplyr::filter(!driver_subclone_identified)%>%
                          pull(comparator)%>%as.numeric()
  )
  
  if(return =="pval"){
    return(wc_test$p.value) 
  }
  if(return == "stat"){
    return(wc_test$statistic)
  }
})
}

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "****",
           ifelse(p < 0.01,  "***",
                  ifelse(p < 0.05,  "**",
                         ifelse(p < 0.1,   "*", "")))))
}


top_annotation = HeatmapAnnotation(Wilcox_stat = check_stats("stat"),
                                   p_value = anno_text(
                                     p_to_stars(check_stats("pval")),
                                     gp = gpar(fontsize = 10, fontface = "bold"),
                                     rot = 0,
                                     just = "bottom",
                                     location = -1.5
                                   ))


plotting_mat= apply(plotting_mat, 2, as.numeric)
plotting_mat = scale(plotting_mat)

colnames(plotting_mat) = case_when(colnames(plotting_mat) == "Invasive.tumour.size..mm." ~ "Inv_Size",
                                   colnames(plotting_mat) == "Ki67_index" ~ "Ki67_FACS",
                                   colnames(plotting_mat) == "min_ki67" ~ "Ki67_IHC_min",
                                   colnames(plotting_mat) == "max_ki67" ~ "Ki67_IHC_max",
                                   colnames(plotting_mat) == "mean_ki67" ~ "Ki67_IHC_mean",
                                   colnames(plotting_mat) == "subclonal_mut_count" ~ "subclonal_mut_count",
                                   colnames(plotting_mat) == "Tumour_Weight" ~ "Tumour_Weight",
                                   colnames(plotting_mat) == "max_length_neoaj_Rx" ~ "Length_NeoadjRx",
                                   colnames(plotting_mat) == "Pathology.Grading.Components_M" ~ "Mitotic_Rate",
                                   colnames(plotting_mat) == "TNM8_Stage" ~ "Stage",
                                   colnames(plotting_mat) == "range_ki67" ~ "Ki67_IHC_range",
                                   colnames(plotting_mat) == "driver_count" ~ "driver_count",
                                   colnames(plotting_mat) == "total_count" ~ "total_count")

max_stat = max(check_stats("stat"))
top_anno = HeatmapAnnotation(
  Wilcox_Test = anno_simple(check_stats("stat"), 
                            pch = p_to_stars(check_stats("pval")),
                            gp  = gpar(col = "white"),
                            pt_gp  = gpar(col = "white"),
                            circlize::colorRamp2(c(0, max_stat), c("#f0f0f0", "#00796b"))
  ),
  border = T,
  col = list(Wilcox_Test = circlize::colorRamp2(c(0, max_stat), c("#f0f0f0", "#00796b")))
)

rownames(plotting_mat) = clinical_data$Trial.ID

max(plotting_mat, na.rm = T)
col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("#4575B4", "#f0f0f0", "#D73027"))

rowsplit = if_else(clinical_data$driver_subclone_identified, "Subclonal \n Ki67 Driver", "No Additional \nDriver")


# WRITE PLOT --------------------------------------------------------------
pdf(file.path(OUT_DIR, "8b_ki67_features.pdf"),height = 13, width = 14)
Heatmap(plotting_mat,
        na_col = "white",
        border = "black",
        row_title_gp = gpar(fontface = "bold", fontsize = 15),
        rect_gp = gpar(col = "white", lwd = 1),
        border_gp = gpar(col = "white"),
        name = "Scaled Value",
        top_annotation = top_anno,
        col = col_fun,
        row_split = rowsplit
)
dev.off()


