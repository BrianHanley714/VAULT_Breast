rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(brms)
library(ggnewscale)
library(ggpubr)
library(ComplexHeatmap)
library(cowplot)
library(qs)
library(tidyverse)
library(openxlsx)
library(mutSignatures)
library(gggenes)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")
MATCHED_PT_CHAR_MBTCGA = file.path(BASE, "data","metadata", "matched_patients_characteristics.txt")
IDMAP = file.path(BASE, "data","metadata", "tumouridmap_MB.txt")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
TCGA_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_TCGA_MB.txt")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
CN_SCORES_PATH = file.path(BASE, "data", "metadata", "copy_number_scores.rdata")
PIK3CA_DOM = file.path(BASE, "data", "metadata", "PIK3CAdomains_P42336_EBI_10022025.tsv")

# READ DATA ---------------------------------------------------------------
TCGA_calls = read.delim(TCGA_VARIANTS_PATH)
mbtcga = read.delim(TCGA_VARIANTS_PATH)
clinical_data = read.delim(VAULT_CLIN_PATH)
maf_out = qread(VAULT_VARIANTS_PATH)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]
matched_patients = c(matched_patients_char$PATIENT_ID, tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)])
matched_patients = matched_patients[!is.na(matched_patients)]
data_cnv_scores = qread(CN_SCORES_PATH)
pik3ca_domains = read.delim(PIK3CA_DOM)



cases_w_HEs = unique(maf_out$Trial_ID[maf_out$sample_type == "FFPE" & maf_out$SingReg_Prim_v_Met %in% c("PRIMARY", "RECURRENCE_RS")])
TCGA_calls = TCGA_calls%>%dplyr::filter(study == "TCGA")
maf_out = maf_out%>%dplyr::filter(Tumor_Sample_Barcode != "HF039_FT003")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "annotate_driver_summary.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "plotting_features.R"))

get_ccf_matrix = function(sample_types, pull_data){
  
  ccf_ffpes = lapply(1:length(driver_calls), function(fn){
    driver_calls[[fn]]%>%dplyr::filter(sample_type == sample_types)%>%
      group_by(Hugo_Symbol, vcf_pos, Consensus_annotated, actionable_summary)%>%
      reframe(mean_ccf = mean(t_VAF))%>%
      arrange(-mean_ccf)%>%
      pull(pull_data)
  })
  
  if(pull_data == "Hugo_Symbol" & sample_types == "FFPE"){
    ccf_ffpes = lapply(1:length(driver_calls), function(fn){
      num_ffpes = num_drivers%>%dplyr::filter(sample_type == "FFPE", grepl(trial_ids[fn], Tumor_Sample_Barcode))%>%nrow()
      
      return(
        paste(driver_calls[[fn]]%>%dplyr::filter(sample_type == sample_types)%>%
                group_by(Hugo_Symbol, vcf_pos, Consensus_annotated, actionable_summary)%>%
                reframe(mean_ccf = mean(t_VAF))%>%
                arrange(-mean_ccf)%>%
                pull(pull_data),
              driver_calls[[fn]]%>%dplyr::filter(sample_type == sample_types)%>%
                mutate(num_ffpe = num_ffpes)%>%
                group_by(Hugo_Symbol, vcf_pos, Consensus_annotated, num_ffpe)%>%
                reframe(count = n(), mean_ccf = mean(t_VAF))%>%
                arrange(-mean_ccf)%>%
                mutate(num_samples = paste0("(", count, "/", num_ffpe, ")"))%>%
                pull(num_samples))
      )
    })
  }
  
  max_length = max(sapply(ccf_ffpes, length))
  
  matrix = matrix(nrow = length(ccf_ffpes), ncol = max_length)
  
  for(fn in 1:length(ccf_ffpes)){
    tryCatch(
      {
        matrix[fn,1:length(ccf_ffpes[[fn]])] = ccf_ffpes[[fn]]    
      },
      error = function(e) {
        message("Error on file ", fn, ": ", e$message)
        return(NULL) 
      }
    )
    
  }
  rownames(matrix) = trial_ids
  return(matrix)
}





# Options -----------------------------------------------------------------
basesize = 20
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"
rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"
dark_grey  = "#4D4D4D"
light_grey = "#D9D9D9"

single_cases = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                Consensus_annotated != "Non_Driver")%>%
  dplyr::select(Tumor_Sample_Barcode, t_VAF, Hugo_Symbol)%>%
  dplyr::filter(!duplicated(paste(Tumor_Sample_Barcode, Hugo_Symbol)))%>%
  dplyr::filter(!Hugo_Symbol %in% Hugo_Symbol[duplicated(Hugo_Symbol)])%>%pull(Hugo_Symbol)



plotting_mat = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                Consensus_annotated != "Non_Driver")%>%
  dplyr::select(Tumor_Sample_Barcode, t_VAF, Hugo_Symbol)%>%
  dplyr::filter(!duplicated(paste(Tumor_Sample_Barcode, Hugo_Symbol)))%>%
  mutate(Hugo_Symbol = if_else(Hugo_Symbol %in% single_cases, "further driver", Hugo_Symbol))%>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
  mutate(t_VAF = mean(t_VAF))%>%
  dplyr::filter(!duplicated(paste(Tumor_Sample_Barcode, Hugo_Symbol)))%>%
  pivot_wider(values_from = t_VAF, names_from = Tumor_Sample_Barcode)%>%as.data.frame()

count_mat = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                Consensus_annotated != "Non_Driver")%>%
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)%>%
  dplyr::filter(!duplicated(paste(Tumor_Sample_Barcode, Hugo_Symbol)))%>%
  mutate(Hugo_Symbol = if_else(Hugo_Symbol %in% single_cases, "further driver", Hugo_Symbol))%>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
  mutate(count = n())%>%
  dplyr::filter(!duplicated(paste(Tumor_Sample_Barcode, Hugo_Symbol)))%>%
  pivot_wider(values_from = count, names_from = Tumor_Sample_Barcode)%>%as.data.frame()


rownames(count_mat) = plotting_mat$Hugo_Symbol
count_mat = count_mat%>%dplyr::select(-Hugo_Symbol)%>%as.matrix()
rownames(plotting_mat) = plotting_mat$Hugo_Symbol

plotting_mat = plotting_mat%>%dplyr::select(-Hugo_Symbol)%>%as.matrix()

colnames(plotting_mat) = substr(colnames(plotting_mat), 1, 5)
cases_wo_drivers = substr(ck_included_cases_vault, 1, 5)[!substr(ck_included_cases_vault, 1, 5) %in% colnames(plotting_mat)]
missing_drivers_mat = matrix(ncol = length(cases_wo_drivers), nrow = nrow(plotting_mat))
colnames(missing_drivers_mat) = cases_wo_drivers

plotting_mat = cbind(plotting_mat, missing_drivers_mat)
count_mat = cbind(count_mat, missing_drivers_mat)


oncoprint_mat = plotting_mat
oncoprint_mat[plotting_mat >0] = "yes"
oncoprint_mat[plotting_mat ==0] = NA

row_order = ComplexHeatmap::row_order(ComplexHeatmap::oncoPrint(oncoprint_mat))
column_order = ComplexHeatmap::column_order(ComplexHeatmap::oncoPrint(oncoprint_mat, ))

clinical_data$receptor_status = clinical_data$Receptor.Subtype
clinical_data_hm = clinical_data[match(colnames(plotting_mat), clinical_data$Trial.ID),]
nice_cols <- c("#616161","#4285f4","#db4437","#f4b400","#0f9d58","#ab47bc","#00acc1",
               "#ff7043","#9e9d24","#5c6bc0","#f06292","#00796b","#c2185b","#7e57c2",
               "#03a9f4","#8bc34a","#fdd835","#fb8c00","#8d6e63","#9e9e9e","#607d8b")

subclonal_tmb = maf_out%>%dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault, clonality == "SUBCLONAL")%>%group_by(Tumor_Sample_Barcode)%>%dplyr::count()%>%mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))
subclonal_tmb_hm = subclonal_tmb$n[match(colnames(plotting_mat), subclonal_tmb$Tumor_Sample_Barcode)]

clonal_tmb = maf_out%>%dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault, clonality == "CLONAL")%>%group_by(Tumor_Sample_Barcode)%>%dplyr::count()%>%mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 5))
clonal_tmb_hm = clonal_tmb$n[match(colnames(plotting_mat), clonal_tmb$Tumor_Sample_Barcode)]

mat_tmb_hm = cbind(subclonal_tmb_hm, clonal_tmb_hm)

mat_tmb_hm = log10(mat_tmb_hm)
mat_tmb_hm[is.na(mat_tmb_hm)] = 0




annotate_driver_summary_tcga = function(dataframe){dataframe%>%
    mutate(ONCOGENIC = if_else(grepl("Onco", ONCOGENIC_ONCOKB)|grepl("annotated", CGI.Oncogenic.Summary), ONCOGENIC_ONCOKB, "Non_Driver"))%>%
    mutate(Biomarker_oncoKB = if_else(HIGHEST_LEVEL_ONCOKB!= "" &!is.na(HIGHEST_LEVEL_ONCOKB), "Biomarker", ONCOGENIC),
           Biomarker_oncoKB = if_else(grepl("Non_Driver|Biomarker", Biomarker_oncoKB), Biomarker_oncoKB, "Driver_Other"),
    )%>%
    mutate(ONCOGENICITY = if_else(grepl("Onco", ONCOGENICITY), ONCOGENICITY, "Non_Driver"),
           Biomarker_ACMG = if_else(grepl("1|2", ACTIONABILITY_TIER), "Biomarker", ONCOGENICITY),
           Biomarker_ACMG = if_else(grepl("Non_Driver|Biomarker", Biomarker_ACMG), Biomarker_ACMG, "Driver_Other")
    )%>%
    mutate(Consensus_annotated = if_else(Biomarker_ACMG == "Non_Driver" & Biomarker_oncoKB == "Non_Driver", "Non_Driver",
                                         if_else(Biomarker_ACMG== "Biomarker" & Biomarker_oncoKB != "Biomarker", Biomarker_ACMG,
                                                 if_else(Biomarker_oncoKB== "Biomarker" & Biomarker_ACMG != "Biomarker", Biomarker_oncoKB,
                                                         if_else(Biomarker_oncoKB== "Biomarker" & Biomarker_ACMG == "Biomarker", "Biomarker", "Driver_Other"))))
           
    )%>%
    mutate(
           rec_type = if_else(study == "TCGA", matched_patients_char$molsubtype[match(Tumor_Sample_Barcode, matched_patients_char$PATIENT_ID)], NA))%>%
    mutate(trial_specific = if_else(!is.na(MCG_n_trials_all_bc)&MCG_n_trials_all_bc>0|
                                      !is.na(MCG_n_trial_TN)&MCG_n_trial_TN>0 &rec_type == "TN"|
                                      !is.na(MCG_n_trial_HER2pos)&MCG_n_trial_HER2pos>0 &rec_type == "HER2-overexpressing"|
                                      !is.na(MCG_n_trials_ERpos)&MCG_n_trials_ERpos>0 &rec_type == "Luminal", 
                                    "trial_inclusion_criterion",
                                    "non_trial_inclusion"))%>%
    mutate(
      actionable_summary = if_else(grepl("Biomarker|Driver_Other", Consensus_annotated) & trial_specific == "trial_inclusion_criterion", paste(Consensus_annotated, trial_specific, sep = ";"), Consensus_annotated))}

TCGA_anno = TCGA_calls%>%
  annotate_driver_summary_tcga()%>%
  dplyr::filter(Consensus_annotated != "Non_Driver")%>%
  mutate(changes = case_when(Hugo_Symbol=="MLL"~ "KMT2A",
                             Hugo_Symbol=="MLL2"~ "KMT2B",
                             Hugo_Symbol=="MLL3"~ "KMT2C",
                             Hugo_Symbol== "MLL4"~ "KMT2D",
  ),
  Hugo_Symbol = if_else(is.na(changes), Hugo_Symbol, changes))%>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = unique(TCGA_calls$Tumor_Sample_Barcode)),
         Hugo_Symbol = if_else(!Hugo_Symbol %in% rownames(plotting_mat), "further driver", Hugo_Symbol),
         Hugo_Symbol = factor(Hugo_Symbol, levels = unique(Hugo_Symbol)))%>%
  dplyr::filter(!is.na(Hugo_Symbol))%>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
  dplyr::count(.drop = F)%>%
  group_by(Hugo_Symbol)%>%
  dplyr::count()%>%
  mutate(Frequency = n/length(unique(TCGA_calls$Tumor_Sample_Barcode)))


TCGA_VAF = TCGA_calls%>%
  annotate_driver_summary_tcga()%>%
  dplyr::filter(Consensus_annotated != "Non_Driver")%>%
  mutate(changes = case_when(Hugo_Symbol=="MLL"~ "KMT2A",
                             Hugo_Symbol=="MLL2"~ "KMT2B",
                             Hugo_Symbol=="MLL3"~ "KMT2C",
                             Hugo_Symbol== "MLL4"~ "KMT2D",
  ),
  Hugo_Symbol = if_else(is.na(changes), Hugo_Symbol, changes))%>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = unique(TCGA_calls$Tumor_Sample_Barcode)),
         Hugo_Symbol = if_else(!Hugo_Symbol %in% rownames(plotting_mat), "further driver", Hugo_Symbol),
         Hugo_Symbol = factor(Hugo_Symbol, levels = unique(Hugo_Symbol)),
         VAF = t_alt_count/(t_alt_count + t_ref_count))%>%
  dplyr::filter(!is.na(Hugo_Symbol))%>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
  reframe(mean_VAF = mean(VAF, na.rm = T))%>%
  group_by(Hugo_Symbol)%>%
  reframe(mean_VAF = mean(mean_VAF, na.rm = T))
col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#1c9099", "#f0f0f0", "#762a83"))
left_anno = rowAnnotation(TCGA_VAF = TCGA_VAF$mean_VAF[match(rownames(plotting_mat), TCGA_VAF$Hugo_Symbol)],
                          TCGA = anno_barplot(TCGA_anno$Frequency[match(rownames(plotting_mat), TCGA_anno$Hugo_Symbol)]),
                          col = list(TCGA_VAF = col_fun))


purity_hm = data_cnv_scores$purity[match(colnames(plotting_mat), data_cnv_scores$Patient_ID)]
ploidy_hm = data_cnv_scores$ploidy[match(colnames(plotting_mat), data_cnv_scores$Patient_ID)]
top_anno = HeatmapAnnotation(Log10_MutCount = anno_barplot(mat_tmb_hm, gp = gpar(border = "white", fill = 1:2, col = "black")),
                             Genomic_Purity = purity_hm,
                             Genomic_Ploidy = ploidy_hm,
                             # Receptor_Subtype = clinical_data_hm$receptor_status, 
                             Histo_Diagnosis = clinical_data_hm$histodx,
                             NeoadjRx = !is.na(clinical_data_hm$neo_adj_treatment_class), 
                             col = list(Receptor_Subtype = c("ER+/HER2-" = "#7e57c2", "HER2+" = "#00796b", "TN" = "#fdd835"),
                                        NeoadjRx = c("TRUE" = nice_cols[1], "FALSE" = "white"),
                                        Histo_Diagnosis = c("IDC_NST" = nice_cols[2], "ILC" = nice_cols[3], "MixedILC/IDC" = nice_cols[4], "MucinousCarcinoma" = nice_cols[5]),
                                        Genomic_Purity = circlize::colorRamp2(c(0, max(purity_hm, na.rm = T)), c("#f0f0f0", nice_cols[5])),
                                        Genomic_Ploidy = circlize::colorRamp2(c(0, max(ploidy_hm, na.rm = T)), c("#f0f0f0", nice_cols[6]))
                             ))

count_mat[is.na(count_mat)] = ""
count_mat[count_mat == "1"] = ""


row_order = row_order[row_order != which(rownames(plotting_mat) == "further driver")]
row_order = c(row_order, which(rownames(plotting_mat) == "further driver"))

rownames(plotting_mat)[rownames(plotting_mat) == "further driver"] = "RARE DRIVERS"

row_split = if_else(rownames(plotting_mat) == "RARE DRIVERS", "rare_genes", "common_genes")

hm_cohort_plot = Heatmap(plotting_mat, col = col_fun,  
                         row_order = row_order, 
                         column_order = column_order, name = "VAF", na_col = "white",border = "black", rect_gp = gpar(colour = "white", border = "grey"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(
            sprintf("%s", count_mat[i, j], digits = 0),
            x, y,
            gp = gpar(fontsize = 10)
          )
        },
        row_split = row_split,
        
        column_split = clinical_data_hm$receptor_status, heatmap_legend_param = list(legend_height = unit(10, "cm")), column_title_gp = gpar(fontface = "bold"),
        top_annotation = top_anno, row_names_side = "left",
        right_annotation = left_anno)


clinical_data$subtype = case_when(
  clinical_data$Receptor.Subtype == "ER+/HER2-" ~ "Luminal",
  clinical_data$Receptor.Subtype == "HER2+" ~ "HER2-overexpressing",
  clinical_data$Receptor.Subtype == "TN" ~ "TN")



# Get driver clonality ----------------------------------------------------


combined_df = 
  bind_rows(TCGA_calls%>%mutate(t_VAF = t_alt_count/(t_alt_count+t_ref_count))%>%
              annotate_driver_summary()%>%dplyr::select(study, Tumor_Sample_Barcode, t_VAF,  Consensus_annotated, clonality, ccf_expected_copies),
            maf_out%>%dplyr::filter( !grepl("common_variant", FILTER.x), 
                                     !Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "Yes"]),
                                     !grepl("HF079", Tumor_Sample_Barcode),
                                     Tumor_Sample_Barcode %in% ck_included_cases_vault)%>% dplyr::select(study, Tumor_Sample_Barcode, Consensus_annotated, clonality, ccf_expected_copies, t_VAF))

plot_fx = function(input){
  input+
    geom_violin(aes(fill = study))+
    geom_boxplot(width = 0.3)+
    geom_jitter(width = 0.05, size = 1)+
    theme_classic(base_size = basesize)+
    stat_compare_means(label.x.npc = 0.5, label = "p.signif", label.y.npc = 0.9)+
    scale_fill_manual(values = c("TCGA" = SR_col, "VAULT" = RS_col))+
    #scale_y_log10()+
    theme(legend.position = "none", axis.title.x = element_blank())
}

grid_tcga_v_VAULT = plot_grid(
  combined_df%>%
    dplyr::filter(!grepl("Non_Driver", Consensus_annotated),
    )%>%
    mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = unique(combined_df$Tumor_Sample_Barcode)),
           study = factor(study, levels = unique(combined_df$study)))%>%
    dplyr::count(Tumor_Sample_Barcode,study, .drop = F)%>%
    dplyr::filter(!(grepl("TCGA", Tumor_Sample_Barcode) & study == "VAULT"),
                  !(grepl("^HF", Tumor_Sample_Barcode) & study == "TCGA"))%>%
    ggplot(aes(study, n))%>%plot_fx()+
    scale_y_log10()+
    ggtitle("Driver Variants"),
  
  combined_df%>%
    group_by(Tumor_Sample_Barcode,study)%>%
    reframe(count = n())%>%
    ggplot(aes(study, count))%>%plot_fx()+
    ggtitle("All Variants")+
    scale_y_log10(),
  
  combined_df%>%
    group_by(Tumor_Sample_Barcode,study)%>%
    reframe(VAF = mean(t_VAF, na.rm = T))%>%
    ggplot(aes(study, VAF))%>%plot_fx()+
    ggtitle("All Variants"),
  
  combined_df%>%
    group_by(Tumor_Sample_Barcode,study)%>%
    reframe(ccf = mean(ccf_expected_copies, na.rm = T))%>%
    ggplot(aes(study, ccf))%>%plot_fx()+
    ggtitle("All Variants"),
  
  combined_df%>%
    dplyr::filter(Consensus_annotated != "Non_Driver",
                  grepl("^CL", clonality))%>%
    group_by(Tumor_Sample_Barcode,study)%>%
    reframe(ccf = mean(ccf_expected_copies, na.rm = T))%>%
    ggplot(aes(study, ccf))%>%plot_fx()+
    ggtitle("Clonal Drivers"),
  combined_df%>%
    dplyr::filter(Consensus_annotated != "Non_Driver",
                  grepl("SUBCL", clonality))%>%
    group_by(Tumor_Sample_Barcode,study)%>%
    reframe(ccf = mean(ccf_expected_copies, na.rm = T))%>%
    ggplot(aes(study, ccf))%>%plot_fx()+
    ggtitle("Subclonal Drivers"), labels = LETTERS[1:6],
  nrow = 3
  
)


progressive_disease = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in%ck_included_cases_vault,
                Consensus_annotated != "Non_Driver",
                
  )%>%
  mutate(progression = clinical_data$disease_progression_after_VAULT_sample[match(substr(Tumor_Sample_Barcode, 1, 5), clinical_data$Trial.ID)],
         time_to_recurrence = clinical_data$days_to_first_recurrence[match(substr(Tumor_Sample_Barcode, 1, 5), clinical_data$Trial.ID)]
  )%>%
  dplyr::filter(progression == "YES")%>%
  dplyr::select(Tumor_Sample_Barcode, time_to_recurrence, Trial_ID, Hugo_Symbol,HGVSp_short, actionable_summary, Consensus_annotated, Hugo_Symbol, t_VAF, clonality)%>%
  arrange(Tumor_Sample_Barcode)%>%
  mutate(Hugo_Symbol = paste(Hugo_Symbol, HGVSp_short, sep = "_"))

row_split_levels = progressive_disease%>%dplyr::filter(!duplicated(Trial_ID))%>%arrange(time_to_recurrence)%>%pull(Trial_ID)
progressive_disease$clonality[is.na(progressive_disease$clonality)] = "INDETERMINATE"
left_anno = rowAnnotation(clonality = progressive_disease$clonality, 
                          col = list(
                            clonality = c("CLONAL" = nice_cols[5],
                                          "SUBCLONAL" = nice_cols[6],
                                          "INDETERMINATE"= nice_cols[7])
                          ),
                          width = unit(ncol(plotting_mat)*10, "mm"),
                          annotation_name_rot = 45,
                          border = T,annotation_name_side = "top", 
                          #show_legend = F,
                          gp = gpar(col = "white", lwd = 1) )
right_anno = rowAnnotation(Actionability = progressive_disease$actionable_summary,
                           #Target_Adj = rep("FALSE", nrow(progressive_disease)),
                           Days_to_recurrence = anno_barplot(progressive_disease$time_to_recurrence),
                           col = list(
                             Actionability = c(
                               "Driver_Other;trial_inclusion_criterion" = nice_cols[1],
                               "Biomarker;trial_inclusion_criterion" = nice_cols[2],   
                               "Driver_Other" = nice_cols[3],
                               "Biomarker" = nice_cols[4]
                             ),
                             Target_Adj = c("FALSE" = "grey") 
                           ),border = T,
                           show_legend = c(Actionability = TRUE, Days_to_recurrence = T,Target_Adj = FALSE),
                           annotation_name_side = "top",
                           annotation_name_rot = 45,
                           gp = gpar(col = "white", lwd = 1) )
plotting_mat = matrix(progressive_disease$t_VAF)
rownames(plotting_mat) = progressive_disease$Hugo_Symbol
colnames(plotting_mat) = "VAF"
col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#1c9099", "#f0f0f0", "#762a83"))


hm = Heatmap(plotting_mat,cluster_rows = F,
             column_names_side = "top",
             col = col_fun,
             heatmap_legend_param = list(
               direction = "vertical"
             ),
             show_row_dend = F, name = "VAF", border = "black",
             rect_gp = gpar(col = "white"),
             row_title_gp = gpar(fontface = "bold"),
             row_names_gp = gpar(fontsize = 6),
             row_split = factor(progressive_disease$Trial_ID, levels = row_split_levels),
             width = unit(ncol(plotting_mat)*6, "mm"),
             left_annotation = left_anno,row_names_side = "left",
             right_annotation = right_anno,column_names_rot = 45,
             row_title_rot = 0)


hm_progressive_disease = ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")

ck_primaryFFPE = maf_out%>%
  dplyr::filter(Trial_ID %in% cases_w_HEs,
                if_else(sample_type == "ENRICHED", Tumor_Sample_Barcode %in% ck_included_cases_vault, TRUE),
                if_else(sample_type == "FFPE", SingReg_Prim_v_Met %in% c("PRIMARY", "RECURRENCE_RS"), TRUE),
                !grepl("_RT|UC", Tumor_Sample_Barcode))

get_num_drivers = TRUE
if(get_num_drivers == TRUE){
  num_drivers = ck_primaryFFPE%>%
    dplyr::filter(!grepl("Non_Driver", Consensus_annotated))%>%
    group_by(Tumor_Sample_Barcode)%>%
    dplyr::count()
  
  num_drivers = rbind(num_drivers,
                      bind_cols(Tumor_Sample_Barcode = ck_primaryFFPE%>%
                                  dplyr::filter(!(Tumor_Sample_Barcode %in% c(num_drivers$Tumor_Sample_Barcode)))%>%pull(Tumor_Sample_Barcode)%>%unique(), n = 0)
  )
  num_drivers = num_drivers%>%mutate(Trial_ID = substr(Tumor_Sample_Barcode, 1, 5),
                                     sample_type = if_else(grepl("ET", Tumor_Sample_Barcode), "ENRICHED", "FFPE"))
  num_drivers
}

trial_ids = unique(ck_primaryFFPE$Trial_ID)
driver_calls = lapply(trial_ids, function(fn){
  return(ck_primaryFFPE%>%
           dplyr::filter(!grepl("Non_Driver", Consensus_annotated))%>%
           dplyr::filter(Trial_ID == fn)%>%
           dplyr::select(Tumor_Sample_Barcode,  t_VAF,vcf_pos, Hugo_Symbol, sample_type, ccf_expected_copies, Consensus_annotated, actionable_summary))
})

plotheatmap = T
if(plotheatmap == T){
  mat_FFPES = get_ccf_matrix(sample_type = "FFPE", pull_data = "mean_ccf")
  mat_ENRICHED = get_ccf_matrix(sample_type = "ENRICHED", pull_data = "mean_ccf")
  
  exclude = rownames(mat_ENRICHED)[rowSums(is.na(mat_ENRICHED)) == ncol(mat_ENRICHED)] # checked these. HF79 has not RepSeq solution, HF390 and HF032 have no
  mat_ENRICHED = mat_ENRICHED[!(rownames(mat_ENRICHED) %in% exclude),]
  mat_FFPES = mat_FFPES[!(rownames(mat_FFPES) %in% exclude),]
  
  HugoSymbols_FFPES = get_ccf_matrix(sample_type = "FFPE", pull_data = "Hugo_Symbol")
  HugoSymbols_ENRICHED = get_ccf_matrix(sample_type = "ENRICHED", pull_data = "Hugo_Symbol")
  HugoSymbols_ENRICHED = HugoSymbols_ENRICHED[!(rownames(HugoSymbols_ENRICHED) %in% exclude),]
  HugoSymbols_FFPES = HugoSymbols_FFPES[!(rownames(HugoSymbols_FFPES) %in% exclude),]
  
  
  
  Consensus_annotated_FFPES = get_ccf_matrix(sample_type = "FFPE", pull_data = "actionable_summary")
  Consensus_annotated_ENRICHED = get_ccf_matrix(sample_type = "ENRICHED", pull_data = "actionable_summary")
  Consensus_annotated_ENRICHED = Consensus_annotated_ENRICHED[!(rownames(Consensus_annotated_ENRICHED) %in% exclude),]
  Consensus_annotated_FFPES = Consensus_annotated_FFPES[!(rownames(Consensus_annotated_FFPES) %in% exclude),]
  
  roworder = names(rowSums(is.na(mat_ENRICHED))%>%sort)
  order = roworder
  roworder = c(paste0("CKsort_", roworder), paste0("FFPE_", roworder))
  
  
  rownames(HugoSymbols_FFPES) = paste0("FFPE_", rownames(HugoSymbols_FFPES))
  rownames(HugoSymbols_ENRICHED) = paste0("CKsort_", rownames(HugoSymbols_ENRICHED))
  
  rownames(mat_FFPES) = paste0("FFPE_", rownames(mat_FFPES))
  rownames(mat_ENRICHED) = paste0("CKsort_", rownames(mat_ENRICHED))
  
  if(ncol(mat_ENRICHED) != ncol(mat_FFPES)){
    difference = ncol(mat_ENRICHED) - ncol(mat_FFPES)
    abs_dff = abs(difference)
    if(difference<0){
      mat_ENRICHED = cbind(mat_ENRICHED, matrix(NA, ncol = abs_dff, nrow = nrow(mat_ENRICHED)))
      Consensus_annotated_ENRICHED= cbind(Consensus_annotated_ENRICHED, matrix(NA, ncol = abs_dff, nrow = nrow(Consensus_annotated_ENRICHED)))
      HugoSymbols_ENRICHED= cbind(HugoSymbols_ENRICHED, matrix(NA, ncol = abs_dff, nrow = nrow(HugoSymbols_ENRICHED)))
    }
    
    if(difference>0){
      mat_FFPES = cbind(mat_FFPES, matrix(NA, ncol = abs_dff, nrow = nrow(mat_FFPES)))
      Consensus_annotated_FFPES= cbind(Consensus_annotated_FFPES, matrix(NA, ncol = abs_dff, nrow = nrow(Consensus_annotated_FFPES)))
      HugoSymbols_FFPES= cbind(HugoSymbols_FFPES, matrix(NA, ncol = abs_dff, nrow = nrow(HugoSymbols_FFPES)))
    }
  }
  
  mat_out = rbind(mat_ENRICHED, mat_FFPES)
  Consensus_annotated_out = rbind(Consensus_annotated_ENRICHED, Consensus_annotated_FFPES)
  HugoSymbols_out = rbind(HugoSymbols_ENRICHED, HugoSymbols_FFPES)
  
  max_value = 1
  
  Heatmap_annotation = rowAnnotation(`Sample Type` = sub("_.*", "", roworder), annotation_name_gp = gpar(col = NA), 
                                     annotation_legend_param = list(
                                       by_row = TRUE,
                                       nrow = 1),
                                     col = list(`Sample Type` = c("CKsort" = "#6A3D9Adb", "FFPE" = "#33A02Cb7")))
  
  low = 
    mid = "white"
  high = 
    
    
  
  "#D95F02"
  
 
  
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#4575B4", "#f0f0f0", "#D73027"))
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#F1C40F", "#f0f0f0",  "#2C3E50"))
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#4C6A92", "#f0f0f0",  "#E9C46A"))
  colnames(mat_out) = 1:ncol(mat_out)
  
  hm_variants_sr_v_rs = Heatmap(mat_out, cluster_rows = F, cluster_columns = F, 
               col = col_fun,
               show_column_names = F,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val <- HugoSymbols_out[i, j]
                 test_val = mat_out[i, j]
                 anno_val = Consensus_annotated_out[i,j]
                 txt_col <- ifelse(test_val > max_value*0.75 | test_val < max_value*0.25, "white", "black")
                 grid.text(sprintf("%s", val), x, y, gp = gpar(fontsize = 9, col = txt_col))
                 if (!is.na(test_val) && grepl("Biomarker|trial_inclusion_criterion", anno_val)) {   
                   grid.rect(x = x, y = y, 
                             width = width, height = height,
                             gp = gpar(col = "black", fill = NA, lwd = 2))
                 }
               },name = "VAF",
               column_order = rev(colnames(mat_out)),
               row_split = factor(sub(".*_", "", rownames(mat_out)), levels = rev(order)),
               na_col = "white",
               show_row_names = F,
               row_title_gp = gpar(fontface = "bold", fontsize = 15),
               rect_gp = gpar(col = "white", lwd = 1),
               row_title_rot = 0, row_gap = unit(3, 'mm'),
               row_title_side = "right",
               heatmap_legend_param = list(
                 legend_width = unit(10, "cm"),
                 direction = "horizontal",
                 legend_height = unit(10, "cm"),
                 grid_width   = unit(10, "mm")
               ), right_annotation = Heatmap_annotation
  )
  

}

ITH_plot = maf_out%>%
  dplyr::filter(!grepl("comm", FILTER.x))%>%
  dplyr::filter(Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "No"]))%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
  mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
  group_by(Tumor_Sample_Barcode, clonality, neoadj)%>%
  dplyr::count()%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  pivot_wider(values_from = "n", names_from = "clonality")%>%
  mutate(ITH = SUBCLONAL/CLONAL)%>%
  ggplot(aes(neoadj, ITH))+
  geom_violin(aes(fill = neoadj), width = 0.8, trim = F, scale = "width")+
  geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
  geom_jitter(width = 0.1, alpha = 0.6)+
  stat_compare_means(label = "p.signif")+
  theme_classic(base_size = basesize)+
  scale_y_log10()+
  scale_fill_manual(values = c(RS_col, SR_col))+
  theme(legend.position = "none", axis.title.x = element_blank())+
  ylab("Intratumour Heterogeneity Index")

SMB_by_treatment_status = 
plot_grid(
  maf_out%>%
    dplyr::filter(!grepl("comm", FILTER.x))%>%
    dplyr::filter(Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "No"]))%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
    mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
    group_by(Tumor_Sample_Barcode, clonality, neoadj)%>%
    dplyr::count()%>%
    dplyr::filter(grepl("CLON", clonality))%>%
    pivot_wider(values_from = "n", names_from = "clonality")%>%
    ggplot(aes(neoadj, SUBCLONAL, label = Tumor_Sample_Barcode))+
    geom_violin(aes(fill = neoadj), width = 0.8, trim = F, scale = "width")+
    geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
    geom_jitter(width = 0.1, alpha = 0.6)+
    stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
    theme_classic(base_size = basesize)+
    scale_y_log10()+
    scale_fill_manual(values = c(RS_col, RS_col))+
    theme(legend.position = "none", axis.title.x = element_blank())+
    ylab("Subclonal Mutation Count")+
    ggtitle("RepSeq"),
  
  
  maf_out%>%
    dplyr::filter(grepl("_B|FT", Tumor_Sample_Barcode))%>%
    mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
    group_by(Tumor_Sample_Barcode, clonality, neoadj, Trial_ID)%>%
    dplyr::count()%>%
    dplyr::filter(grepl("CLON", clonality))%>%
    pivot_wider(values_from = "n", names_from = "clonality")%>%
    ggplot(aes(neoadj, SUBCLONAL, label = Tumor_Sample_Barcode))+
    geom_violin(aes(fill = neoadj), width = 0.8, trim = F, scale = "width")+
    geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
    geom_jitter(width = 0.1, alpha = 0.6)+
    stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
    theme_classic(base_size = basesize)+
    scale_y_log10()+
    scale_fill_manual(values = c(SR_col, SR_col))+
    theme(legend.position = "none", axis.title.x = element_blank())+
    ylab("Subclonal Mutation Count")+
    ggtitle("FFPE Samples")
)

plot_by_treatment = function(treatment = "Naive", plot_y_axis = "SUBCLONAL", plot_y_axis_label = "PROVIDE Y LABEL"){
  
  maf_out%>%
    dplyr::filter(Variant_Classification != "Silent")%>%
    dplyr::filter(Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "No"]))%>%
    mutate(TStage = clinical_data$T.Stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
    mutate(NStage = clinical_data$N.stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
    mutate(MStage = clinical_data$M.stage[match(Trial_ID, clinical_data$Trial.ID)],
           MStage = if_else(is.na(MStage), "unknown", MStage)
           )%>%
    mutate(Subtype = clinical_data$subtype[match(Trial_ID, clinical_data$Trial.ID)])%>%
    mutate(Inv_tumour_size = clinical_data$Invasive.tumour.size..mm.  [match(Trial_ID, clinical_data$Trial.ID)])%>%
    mutate(histodx = clinical_data$Histological.Subtype  [match(Trial_ID, clinical_data$Trial.ID)])%>%
    mutate(Recurrence = !is.na(clinical_data$days_to_first_recurrence)[match(Trial_ID, clinical_data$Trial.ID)])%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
    mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
    dplyr::filter(neoadj == treatment)%>%
    group_by(Tumor_Sample_Barcode,Inv_tumour_size, clonality, neoadj, TStage, NStage, Recurrence, MStage, histodx)%>%
    dplyr::count()%>%
    dplyr::filter(grepl("CLON", clonality))%>%
    pivot_wider(values_from = "n", names_from = "clonality")%>%
    dplyr::filter(NStage!="Nx")%>%
    mutate(N_low = grepl("N0", NStage),
           test = 
             Recurrence | MStage == "M1"
           | ! N_low
    )%>%
    mutate(ITH = SUBCLONAL/CLONAL)%>%
    ggplot(aes(test, .data[[plot_y_axis]], label = Tumor_Sample_Barcode))+
    geom_violin(aes(fill = test))+
    geom_boxplot(width = .2, outlier.alpha = 0)+
    stat_compare_means(label = "p.signif", label.x.npc = 0.5, label.y.npc = 0.8)+
    scale_y_log10()+
    ylab(plot_y_axis_label)+
    geom_jitter(width = 0.1)+
    theme_classic(base_size = basesize)+
    scale_fill_manual(values = c(RS_col, RS_col))+
    theme(legend.position = "none")+
    xlab("Metastatic Disease")+
    ggtitle(treatment)
  
}


maf_out%>%
  dplyr::filter(Variant_Classification != "Silent")%>%

  dplyr::filter(Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "No"]))%>%
  mutate(TStage = clinical_data$T.Stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(NStage = clinical_data$N.stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(MStage = clinical_data$M.stage[match(Trial_ID, clinical_data$Trial.ID)],
         MStage = if_else(is.na(MStage), "unknown", MStage)
  )%>%
  mutate(Subtype = clinical_data$subtype[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Inv_tumour_size = clinical_data$Invasive.tumour.size..mm.  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(histodx = clinical_data$Histological.Subtype  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Recurrence = !is.na(clinical_data$days_to_first_recurrence)[match(Trial_ID, clinical_data$Trial.ID)])%>%
  dplyr::filter(Tumor_Sample_Barcode %in% c(ck_primaryFFPE$Tumor_Sample_Barcode),
                !Tumor_Sample_Barcode %in% ck_included_cases_vault
                  )%>%
  mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
  dplyr::filter(neoadj == "Naive")%>%
  group_by(Tumor_Sample_Barcode,Inv_tumour_size, clonality, neoadj, TStage, NStage, Recurrence, MStage, histodx)%>%
  dplyr::count()%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  pivot_wider(values_from = "n", names_from = "clonality")%>%
  mutate(N_low = grepl("N0", NStage),
         test = 
           Recurrence | MStage == "M1"
         | ! N_low
  )%>%
  mutate(ITH = SUBCLONAL/CLONAL)%>%
  ggplot(aes(test, SUBCLONAL, label = Tumor_Sample_Barcode))+
  geom_violin(aes(fill = test))+
  geom_boxplot(width = .2, outlier.alpha = 0)+
  stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
  scale_y_log10()+
  geom_jitter(aes(col = Inv_tumour_size), width = 0.1)+
  theme_classic(base_size = basesize)+
  scale_fill_manual(values = c(RS_col, SR_col))+
  theme(legend.position = "none")+
  xlab("Metastatic Disease")+
  ylab("Subclonal Mutation Count")+
  ggtitle("FFPE Samples")

disease_risk_RS = maf_out%>%
  dplyr::filter(Variant_Classification != "Silent")%>%
  mutate(TNM8Stage = clinical_data$TNM8_Stage[match(Trial_ID, clinical_data$Trial.ID)],
  )%>%
  mutate(Subtype = clinical_data$subtype[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Inv_tumour_size = clinical_data$Invasive.tumour.size..mm.  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(histodx = clinical_data$Histological.Subtype  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Recurrence = (clinical_data$recurrence_current_sample_is_recurence[match(Trial_ID, clinical_data$Trial.ID)]) == "Yes")%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault
  )%>%
  mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
  group_by(Tumor_Sample_Barcode,Inv_tumour_size, clonality, neoadj, TNM8Stage, Recurrence, histodx)%>%
  dplyr::count()%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  pivot_wider(values_from = "n", names_from = "clonality")%>%
  mutate(test = 
           if_else(Recurrence | TNM8Stage %in% c(3, 4), "high_risk", "low_risk"),
         test = factor(test, levels = c("low_risk", "high_risk"))
  )%>%
  mutate(ITH = SUBCLONAL/CLONAL)%>%
  ggplot(aes(test, SUBCLONAL))+
  geom_violin(aes(fill = test))+
  geom_boxplot(width = .2, outlier.alpha = 0)+
  stat_compare_means(
    label = "p.signif", label.x.npc = 0.5, label.y.npc = 0.8 
    )+
  scale_y_log10()+
  geom_jitter(width = 0.1)+
  theme_classic(base_size = basesize)+
  scale_fill_manual(values = c(RS_col, RS_col))+
  theme(legend.position = "none")+
  ylab("Subclonal Mutation Count")+
  theme(axis.title.x = element_blank())+
  ggtitle("RepSeq")


disease_risk_SOC = maf_out%>%
  dplyr::filter(Variant_Classification != "Silent")%>%
  mutate(TNM8Stage = clinical_data$TNM8_Stage[match(Trial_ID, clinical_data$Trial.ID)],
  )%>%
  mutate(Subtype = clinical_data$subtype[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Inv_tumour_size = clinical_data$Invasive.tumour.size..mm.  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(histodx = clinical_data$Histological.Subtype  [match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Recurrence = (clinical_data$recurrence_current_sample_is_recurence[match(Trial_ID, clinical_data$Trial.ID)]) == "Yes")%>%
  dplyr::filter(Tumor_Sample_Barcode %in% c(ck_primaryFFPE$Tumor_Sample_Barcode),
               !Tumor_Sample_Barcode %in% ck_included_cases_vault
  )%>%
  mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
  group_by(Tumor_Sample_Barcode,Inv_tumour_size, clonality, neoadj, TNM8Stage, Recurrence, histodx)%>%
  dplyr::count()%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  pivot_wider(values_from = "n", names_from = "clonality")%>%
  mutate(test = 
           if_else(Recurrence | TNM8Stage %in% c(3, 4), "high_risk", "low_risk"),
         test = factor(test, levels = c("low_risk", "high_risk"))
  )%>%dplyr::filter(!is.na(test))%>%
  ggplot(aes(test, SUBCLONAL))+
  geom_violin(aes(fill = test))+
  geom_boxplot(width = .2, outlier.alpha = 0)+
  stat_compare_means(
    label = "p.signif", label.x.npc = 0.5, label.y.npc = 0.8
  )+
  scale_y_log10()+
  geom_jitter(width = 0.1)+
  theme_classic(base_size = basesize)+
  scale_fill_manual(values = c(SR_col, SR_col))+
  theme(legend.position = "none")+
  ylab("Subclonal Mutation Count")+
  theme(axis.title.x = element_blank())+
  ggtitle("FFPE Samples")

actionability_out = bind_rows(
  maf_out%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                  !is.na(HIGHEST_LEVEL))%>%
    mutate(framework = "ONCOKB")%>%
    dplyr::select(ccf_expected_copies, Tier = HIGHEST_LEVEL, framework, clonality),
  
  maf_out%>%
    dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                  !is.na(ACTIONABILITY_TIER))%>%
    mutate(framework = "ACMG",
           ACTIONABILITY_TIER = as.character(ACTIONABILITY_TIER),
           ACTIONABILITY_TIER =  case_when(
             ACTIONABILITY_TIER == "1"~ "Tier I",
             ACTIONABILITY_TIER == "2"~ "Tier II",
             ACTIONABILITY_TIER == "3"~ "Tier III",
             ACTIONABILITY_TIER == "4"~ "Tier IV",
             ACTIONABILITY_TIER == "5"~ "Tier V"
           )
    )%>%
    dplyr::select(ccf_expected_copies, Tier = ACTIONABILITY_TIER, framework, clonality)
)

shade_between_white <- function(col, n = 5) {
  target = grDevices::col2rgb(col) / 255
  
  white = c(1, 1, 1)
  
  shades = t(sapply(seq(0, 1, length.out = n), function(a) {
    (1 - a) * white + a * target
  }))

  grDevices::rgb(shades[,1], shades[,2], shades[,3])
}

blue_shades = shade_between_white("navy", 5)



non_synonymous = c("Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "In_Frame_Ins")  

clinical_trial_biomarkers_ccf_plot = plot_grid(
actionability_out%>%
  dplyr::filter(!Tier %in% c("Tier IV", "Tier V"))%>%
  ggplot(aes(framework, ccf_expected_copies))+
  geom_violin(aes(), width = 0.8, trim = T, scale = "width")+
  geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
  geom_jitter(aes(fill= Tier, col = Tier), width = 0.05, vjust = 0, shape = 21, size = 5)+
  scale_color_manual(name = "OncoKB", values= c("LEVEL_1" = "black", "LEVEL_2" ="black", "LEVEL_3A" = "black", "LEVEL_3B" = "black", "LEVEL_4" = "black"), na.value = "#0000FF00")+
  scale_fill_manual(name = "OncoKB",  values= c("LEVEL_1" = blue_shades[5], "LEVEL_2" = blue_shades[4], "LEVEL_3A" = blue_shades[3], "LEVEL_3B" = blue_shades[2], "LEVEL_4" = blue_shades[1]),
                    na.value = "#0000FF00")+
  new_scale_fill()+
  new_scale_color()+
  geom_jitter(aes(fill= Tier, col = Tier), width = 0.05, vjust = 0, shape = 21, size = 5)+
  scale_color_manual(name = "ACMG", values= c("Tier I" = "black", "Tier II" ="black", "Tier III" = "black"), na.value = "#0000FF00")+
  scale_fill_manual(name = "ACMG", values= c("Tier I" = blue_shades[5], "Tier II" = blue_shades[3], "Tier III" = blue_shades[1]), na.value = "#0000FF00")+
  theme_classic(base_size = basesize)+
  ylab("Cancer Cell Fraction")+
  theme(axis.title.x = element_blank(),
        legend.position = c(0.5, 0.4)),
maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, ccf_expected_copies, MCG_n_trials_all_bc, MCG_n_trials_ERpos, MCG_n_trial_HER2pos, MCG_n_trial_TN, Variant_Classification, CLINVAR_CLASSIFICATION, ONCOGENICITY, ONCOGENIC, rec_type)%>%
  dplyr::filter(ONCOGENICITY != "Non_Driver" | grepl("pathogenic", CLINVAR_CLASSIFICATION) | grepl("Onco", ONCOGENIC) )%>%
  pivot_longer(c(MCG_n_trials_all_bc, MCG_n_trials_ERpos, MCG_n_trial_HER2pos, MCG_n_trial_TN), values_to = "n_trials")%>%
  dplyr::filter(!is.na(n_trials) & n_trials != "0", 
                !(name == "MCG_n_trial_TN" & rec_type != "TN"),
                !(name == "MCG_n_trial_HER2pos" & rec_type != "HER2+"),
                !(name == "MCG_n_trials_ERpos" & rec_type != "ER+/HER2-")
  )%>%
  mutate(n_trials = as.numeric(n_trials),
         name = case_when(name == "MCG_n_trial_TN" ~ "Triple Negative",
                          name == "MCG_n_trial_HER2pos" ~ "HER2+",
                          name == "MCG_n_trials_ERpos" ~ "ER+/HER2-",
                          name == "MCG_n_trials_all_bc" ~ "Subtype Agnostic"
         ),
         name = factor(name, levels = c("Triple Negative", "HER2+", "ER+/HER2-", "Subtype Agnostic"))
  )%>%
  ggplot(aes(name, ccf_expected_copies))+
  geom_violin(aes(), width = 0.8, trim = T, scale = "width")+
  geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
  geom_jitter(aes(fill = n_trials), width = 0.05, shape = 21, size = 5)+
  scale_fill_gradient(low = "white", high = "maroon")+
  theme_classic(base_size = basesize)+
  rot.lab()+
  ylab("Cancer Cell Fraction")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = c(0.2, 0.3))+
  ggtitle("Clinical Trial - Inclusion Criteria"),
nrow = 2)



# Lollipop plot -----------------------------------------------------------
vault = maf_out
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%mutate(MCG_n_trials_all_bc = as.integer(MCG_n_trials_all_bc),
                                       MCG_n_trials_ERpos = as.integer(MCG_n_trials_ERpos),
                                       MCG_n_trial_HER2pos = as.integer(MCG_n_trial_HER2pos),
                                       MCG_n_trial_TN = as.integer(MCG_n_trial_TN)
                        )%>% dplyr::select(common_names)
)
# PIK3CA

pik3ca_domains = pik3ca_domains%>%
  dplyr::filter(Source.Database == "pfam")%>%
  mutate(start = as.numeric(sub("\\.\\..*", "", Matches)),
         end = as.numeric(sub(".*\\.\\.", "", Matches)),
         gene = "PIK3CA")%>%
  mutate(Name = case_when(Name == "Phosphatidylinositol 3- and 4-kinase" ~"Kinase",
                          Name == "Phosphoinositide 3-kinase family, accessory domain (PIK domain)"~"Accessory",
                          Name == "Phosphoinositide 3-kinase C2"~ "C2",
                          Name == "PI3-kinase family, ras-binding domain"~ "ras-binding",
                          Name == "PI3-kinase family, p85-binding domain"~ "p85-binding")
  )
mutations = combined_df%>%
  dplyr::filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%
  dplyr::filter(if_else(study == "VAULT", Tumor_Sample_Barcode %in% ck_included_cases_vault, TRUE))%>%
  mutate(study = if_else(study %in% c("METABRIC", "TCGA"), "SingReg", "RepSamp"))%>%
  dplyr::filter(Hugo_Symbol == "PIK3CA")%>%
  dplyr::select(AMINO_ACID_START, study, clonality, HGVSp_short) 

# PIK3CA lollipop ---------------------------------------------------------


mutations = mutations%>%
  group_by(study, clonality, AMINO_ACID_START)%>%
  reframe(Frequency = n())%>%
  ungroup()%>%
  mutate(AMINO_ACID_START = as.numeric(AMINO_ACID_START))%>%
  mutate(Frequency = if_else(clonality == "SUBCLONAL", Frequency - (2*Frequency), Frequency),
         study = factor(study, levels = c("SingReg", "RepSamp")),
         clonality =if_else(study == "SingReg", "TCGA/METABRIC", clonality),
         clonality = factor(clonality, levels = c("TCGA/METABRIC", 
                                                  "CLONAL", 
                                                  "SUBCLONAL")),
         Frequency2 = if_else(clonality == "TCGA/METABRIC", NA, Frequency)
  )




# DRAW PLOT ---------------------------------------------------------------
pik3ca_lollipop = ggplot() +
  geom_segment(data = mutations, aes(x = AMINO_ACID_START, xend = AMINO_ACID_START, y = 0, yend = Frequency, color = clonality), size = 1) +
  geom_point(data = mutations, aes(x = AMINO_ACID_START, y = Frequency, color = clonality, fill = clonality), 
             size = 1, shape = 21,  stroke = 1.2) +
  geom_segment(data = mutations, aes(x = AMINO_ACID_START, xend = AMINO_ACID_START, y = 0, yend = Frequency2, color = clonality), size = 1) +
  geom_point(data = mutations, aes(x = AMINO_ACID_START, y = Frequency2, color = clonality, fill = clonality), 
             size = 1, shape = 21,  stroke = 1.2) +
  geom_rect(aes(xmin = 0, xmax = unique(pik3ca_domains$Protein.Length), ymin = -0.1, ymax = 0.1), fill = "grey", col = "black")+
  
  geom_gene_arrow(data = pik3ca_domains, 
                  aes(xmin = start,xmax = end, y = 0, fill = Name),  arrowhead_width = grid::unit(x = 0, "mm"), arrowhead_height = grid::unit(0, "mm")) +
  scale_fill_manual(values = c("Kinase" = colours[1],
                               "Accessory" = colours[2],
                               "C2" = colours[3],
                               "ras-binding" = colours[4],
                               "p85-binding" = colours[5],
                               "CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col,
                               "TCGA/METABRIC" = "lightgrey",
                               "RepSamp" = rs_col),
                    breaks = unique(pik3ca_domains$Name))+
  scale_y_continuous(trans = "pseudo_log", breaks = c(-1, -10, -100, 0, 1, 10, 100)) +
  scale_color_manual(values = c("TCGA/METABRIC" = "lightgrey", "CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col,
                                "SingReg" = "lightgrey",
                                "RepSamp" = rs_col), )+
  theme_classic(base_size = 20)+
  xlab("Amino Acid Position")+
  ylab("Variant Count")+
  ggtitle("PIK3CA")+
  guides(col = "none")+
  labs(fill = NULL)+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  annotate("text", x = unique(pik3ca_domains$Protein.Length)*0.05, y = max(mutations$Frequency), label = "CLONAL", col = clonal_col, fontface = "bold", size = 6)+
  annotate("text", x = unique(pik3ca_domains$Protein.Length)*0.05, y = min(mutations$Frequency), label = "SUBCLONAL", col = subclonal_col, fontface = "bold", size = 6)


# multivariate model  -----------------------------------------------------


test_df = maf_out%>%
  dplyr::filter(!(grepl("Silent", Variant_Classification)))%>%
  dplyr::filter(Trial_ID %in% c(clinical_data$Trial.ID[clinical_data$recurrence_current_sample_is_recurence == "No"]))%>%
  mutate(TStage = clinical_data$T.Stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(NStage = clinical_data$N.stage[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Age = clinical_data$Age[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Histological.grade = clinical_data$Histological.grade[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(LVI = clinical_data$LVI[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Inv_tumour_size = clinical_data$Invasive.tumour.size..mm.[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Histological_Subtype = clinical_data$Histological.Subtype[match(Trial_ID, clinical_data$Trial.ID)])%>%
  mutate(Recurrence = !is.na(clinical_data$days_to_first_recurrence)[match(Trial_ID, clinical_data$Trial.ID)])%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
  mutate(neoadj = if_else(Trial_ID %in% c(clinical_data$Trial.ID[!is.na(clinical_data$neo_adj_treatment_class)]), "Neoadj", "Naive"))%>%
  dplyr::filter(neoadj == "Naive")%>%
  group_by(Tumor_Sample_Barcode,Inv_tumour_size, Age, LVI, Histological_Subtype, Histological.grade, clonality, NStage, neoadj, TStage, Recurrence, rec_type)%>%
  dplyr::count()%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  pivot_wider(values_from = "n", names_from = "clonality")%>%
  mutate(ITH = SUBCLONAL/CLONAL)%>%
  mutate(N_low = grepl("N0", NStage),
         test = Recurrence | ! N_low)



test_df$Histological.grade = if_else(test_df$Histological.grade=="G3", "high", "low")
test_df$Histological_Subtype = if_else(test_df$Histological_Subtype %in% c("IDC_NST", "MucinousCarcinoma"), "IDC", "ILC")
test_df$rec_type = if_else(test_df$rec_type == "ER+/HER2-", "ER+/HER2-", "Non-ER+/HER2-")

fit_ITH <- glm(test ~ Inv_tumour_size+SUBCLONAL+LVI+Age+Histological_Subtype+Histological.grade+SUBCLONAL,
               data = test_df,
               family = binomial
)

summary(fit_ITH)







# multivariate bayesian model ---------------------------------------------


scale_test_df = test_df


scale_test_df[,c("Inv_tumour_size", "SUBCLONAL", "Age")] =  scale(test_df[,c("Inv_tumour_size", "SUBCLONAL", "Age")])

test = brm(test ~ Inv_tumour_size+SUBCLONAL+LVI+Age+Histological_Subtype+Histological.grade+SUBCLONAL,
           data = scale_test_df,
           family = bernoulli(link = "logit"),
           prior = c(
             prior(normal(0, 1), class = "b")

           ),
           chains = 2,
           iter = 4000,
           warmup =1000,
           cores = 1,
           seed = 123
)


post = posterior_samples(test, "b_SUBCLONAL")
p_two_sided <- 2 * min(
  mean(post$b_SUBCLONAL > 0),
  mean(post$b_SUBCLONAL < 0)
)
p_two_sided

hypothesis(test, "Inv_tumour_size > 0")
post = posterior_samples(test, "b_Inv_tumour_size")
p_two_sided <- 2 * min(
  mean(post$b_Inv_tumour_size > 0),
  mean(post$b_Inv_tumour_size < 0)
)
p_two_sided

hypothesis(test, "LVIPresent > 0")
post = posterior_samples(test, "b_LVIPresent")
p_two_sided <- 2 * min(
  mean(post$b_LVIPresent > 0),
  mean(post$b_LVIPresent < 0)
)
p_two_sided
test$data
multivariate_model_out = as.data.frame(test)


descriptive_stats = multivariate_model_out%>%
  pivot_longer(cols = everything())%>%
  group_by(name)%>%
  reframe(Mean = mean(value),
          CI_lower_95 = quantile(value, probs = c(0.025, 0.975))[1],
          CI_upper_95 = quantile(value, probs = c(0.025, 0.975))[2],
          CI_lower_80 = quantile(value, probs = c(0.1, 0.9))[1],
          CI_upper_80 = quantile(value, probs = c(0.1, 0.9))[2]
          
  )%>%dplyr::filter(!name %in% c("lp__", "lprior", "Intercept","b_Intercept"))%>%
  mutate(name = sub("b_", "", name))

unique(descriptive_stats$name)

rename = c("Age" = "Age","Histological.gradelow" = "Low_Grade","Histological.gradeG3"  = "Grade 3", "Histological_SubtypeILC" =  "ILC", "Histological_SubtypeMixedILCDIDC" = "Mixed", 
           "Histological_SubtypeMucinousCarcinoma" = "Mucinous", "Inv_tumour_size" = "Invasive Tumour Size", "LVIPresent"= "Lymphovascular Invasion",
           "SUBCLONAL" = "Subclonal Mutation Burden", "rec_typeERPHER2" = "HER2+","rec_typeTN" = "Triple Negative")                           


levels = c("Subclonal Mutation Burden", "Invasive Tumour Size", "Lymphovascular Invasion","ILC","Low_Grade", "Age")

type = c("Subclonal Mutation Burden" = "Genomic", "Invasive Tumour Size" = "Staging", "Lymphovascular Invasion" = "Histological", "Low_Grade" = "Histological Grade", "Grade 2" = "Histological Grade", "Mixed" =  "Histological Subtype",
         "ILC" =  "Histological Subtype", "Mucinous" =  "Histological Subtype", "Triple Negative" =  "Receptor Subtype", "HER2+" =  "Receptor Subtype", "Age" = "Age")

colours = c("Genomic" = "#8d6e63" ,"Staging"="#f4b400", "Histological" = "#FF7F00db","Histological Grade" = "#1F78B4db" ,"Histological Subtype" = "#E31A1Cff","Receptor Subtype" = "#CAB2D6b7" ,"Age" = "#A6CEE393")

descriptive_stats$name = factor(rename[match(descriptive_stats$name, names(rename))], levels = levels)
descriptive_stats$type = type[match(descriptive_stats$name, names(type))]

linewidth = 2



mv_model = multivariate_model_out%>%
  pivot_longer(cols = everything())%>%
  dplyr::filter(!name %in% c("lp__", "lprior", "Intercept","b_Intercept"))%>%
  
  mutate(name = sub("b_", "", name))%>%
  mutate(name = factor(rename[match(name, names(rename))], levels = levels),
         type = type[match(name, names(type))])%>%
  #mutate(name = if_else(name == "SUBCLONAL", "Subclonal Mutation Burden", name))%>%
  
  ggplot(aes(value))+
  facet_wrap(~name, ncol = 1, scale = "free_y", strip.position = "left")+
  geom_density(aes(fill = type, col = type), alpha = .3)+
  theme_classic(base_size = basesize)+
  theme(
    strip.placement = "outside",  
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  )+
  geom_point(data = descriptive_stats, aes(Mean, 0, col = type), size = linewidth*5)+
  geom_segment(data = descriptive_stats, aes(x = CI_lower_80, xend = CI_upper_80, y = 0, yend = 0, col = type), linewidth = linewidth*3)+
  geom_segment(data = descriptive_stats, aes(col = type, x = CI_lower_95, xend = CI_upper_95, y = 0, yend = 0), linewidth = linewidth)+
  geom_vline(xintercept = 0, linetype = "longdash")+
  #scale_fill_manual(values = colours)+
  scale_fill_manual(values = rep("grey", length(colours)))+
  #scale_colour_manual(values = colours)+
  scale_colour_manual(values = rep(dark_grey, length(colours)))+
  xlab("Standardised Effect Size")+
  guides(fill= guide_legend(nrow = 1, override.aes = list(alpha = 1)), 
         col = "none", 
         linetype = "none",
         shape = "none")



ITH_plot
hm_progressive_disease
ggsave(file.path(OUT_DIR, "3F_recurrent_disease_cohort.png"))
hm_cohort_plot
ggsave(file.path(OUT_DIR, "3A_cohort_overview.png"))
grid_tcga_v_VAULT
ggsave(file.path(OUT_DIR, "SupplementaryFigure_TCGA_v_VAULT.png"))


pdf(file.path(OUT_DIR, "5I_hm_variants_ss.pdf"),height = 13, width = 14)
draw(hm_variants_sr_v_rs, heatmap_legend_side = "bottom")
dev.off()
pik3ca_lollipop
ggsave(file.path(OUT_DIR, "3D_PIK3CA_lollipop.png"))
clinical_trial_biomarkers_ccf_plot
ggsave(file.path(OUT_DIR, "3B_C_clinical_trials_biomarkers_ccf.png"))



plot_grid(
plot_grid(
  plot_by_treatment(treatment = "Naive", plot_y_axis = "CLONAL", plot_y_axis_label = "CLONAL MUTATION COUNT"),
  plot_by_treatment(treatment = "Neoadj", plot_y_axis = "CLONAL", plot_y_axis_label = "CLONAL MUTATION COUNT")
),

plot_grid(
  plot_by_treatment(treatment = "Naive", plot_y_axis = "ITH", plot_y_axis_label = "ITH INDEX"),
  plot_by_treatment(treatment = "Neoadj", plot_y_axis = "ITH", plot_y_axis_label = "ITH INDEX")
),
nrow = 2
)
ggsave(file.path(OUT_DIR, "SupplementaryFigure_clonal_ITH_metastatic_disease.png"))
plot_grid(
plot_grid(
SMB_by_treatment_status,
plot_grid(
disease_risk_RS,
disease_risk_SOC,
nrow = 1),
nrow = 1),
plot_grid(
plot_grid(
  plot_by_treatment(treatment = "Naive", plot_y_axis = "SUBCLONAL", plot_y_axis_label = "Subclonal Mutation Count"),
  plot_by_treatment(treatment = "Neoadj", plot_y_axis = "SUBCLONAL", plot_y_axis_label = "Subclonal Mutation Count")
),
mv_model,
nrow = 1),
nrow = 2
)
ggsave(file.path(OUT_DIR, "4_subclonal_mutation_burden_VAULT.png"))
