rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(survival)
library(survminer)
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
TCGA_PATIENT_CLINICAL_PATH = file.path(BASE, "data", "metadata", "data_clinical_patient.txt")
TCGA_RX_PATH = file.path(BASE, "data", "metadata", "data_timeline_treatment.txt")
DETECTIONS_PATH = file.path(BASE, "data", "image_analysis", "measurements_nuclear_area.tsv")
HF270_EXAMPLE_RS_PATH = file.path(BASE, "data", "variants", "HF270_ET001_D01_facets.RData")
HF270_EXAMPLE_SR_PATH = file.path(BASE, "data", "variants", "HF270_BA013_F01_facets.RData")


# READ DATA ---------------------------------------------------------------
Detections = read.delim(DETECTIONS_PATH)
TCGA_variants = read.delim(TCGA_VARIANTS_PATH)
TCGA_clinical = read.delim(TCGA_PATIENT_CLINICAL_PATH, skip = 4)
TCGA_Rx = read.delim(TCGA_RX_PATH)
clinical_data = read.delim(VAULT_CLIN_PATH)
maf_out = qread(VAULT_VARIANTS_PATH)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]
matched_patients = c(matched_patients_char$PATIENT_ID, tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)])
matched_patients = matched_patients[!is.na(matched_patients)]
data_cnv_scores = qread(CN_SCORES_PATH)
facets = load(HF270_EXAMPLE_RS_PATH)
facets_snps = facets_output$snps
facets_sr = load(HF270_EXAMPLE_SR_PATH)
facets_snps_sr = facets_output$snps


cases_w_HEs = unique(maf_out$Trial_ID[maf_out$sample_type == "FFPE" & maf_out$SingReg_Prim_v_Met %in% c("PRIMARY", "RECURRENCE_RS")])

maf_out = maf_out%>%dplyr::filter(Tumor_Sample_Barcode != "HF039_FT003")


# Graphical Parameters ----------------------------------------------------
basesize=20
RS_col = "#54278F"
SR_col = "#33A02Cb7"
label_size = 6
dark_grey  = "#4D4D4D"
light_grey = "#D9D9D9"

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "annotate_driver_summary.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "plotting_features.R"))


purity_across_cohorts = bind_rows(
  TCGA_variants%>%
    dplyr::filter(Tumor_Sample_Barcode %in% matched_patients_char$PATIENT_ID)%>%
    group_by(Tumor_Sample_Barcode)%>%
    reframe(purity = unique(purity))%>%
    mutate(sample_type = if_else(substr(Tumor_Sample_Barcode, 1, 1) == "T", "TCGA", "METABRIC")),
  maf_out%>%group_by(Tumor_Sample_Barcode, sample_type)%>%
    reframe(purity = unique(purity))%>%
    mutate(sample_type = paste0("VAULT_", sample_type))
)%>%
  dplyr::filter(!grepl("Other", sample_type))%>%
  mutate(sample_type = factor(sample_type, levels = c("VAULT_FFPE", "METABRIC", "TCGA", "VAULT_ENRICHED")))%>%
  ggplot(aes(sample_type, purity))+
  geom_violin(aes(fill = sample_type))+
  geom_boxplot(width = 0.2)+
  geom_point(aes(1, 1.1), alpha = 0)+
  geom_jitter(width = 0.05, alpha = 0.2, size = 1)+
  theme_classic(base_size = basesize)+
  scale_fill_manual(values = c(
    "METABRIC" = light_grey,
    "TCGA" = light_grey,
    "VAULT_FFPE" = SR_col,
    "VAULT_ENRICHED" = RS_col
  ))+
  stat_compare_means(label = "p.signif", label.y = 1, size = label_size, comparisons = list(
    c("METABRIC", "VAULT_FFPE"), c("TCGA", "VAULT_FFPE"),
    c("METABRIC", "VAULT_ENRICHED"), c("TCGA", "VAULT_ENRICHED")
  ))+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
  rot.lab()+
  ylab("Purity")




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


driver_counts_TCGA = TCGA_variants%>%
  annotate_driver_summary_tcga()%>%
  dplyr::filter(
    #!Hugo_Symbol %in% c("PIK3CA", "MAP3K1", "GATA3"),
    Consensus_annotated != "Non_Driver")%>%
  group_by(Tumor_Sample_Barcode, clonality)%>%
  dplyr::count()%>%
  pivot_wider(names_from = clonality, values_from = n)%>%
  mutate(total_driver_count = sum(CLONAL, SUBCLONAL, `NA`, na.rm = T),
         SUBCLONAL = if_else(is.na(SUBCLONAL), 0, SUBCLONAL),
         prop_subclonLdrivers = SUBCLONAL/total_driver_count)

ITH_tmb_TCGA = TCGA_variants%>%
  
  dplyr::filter(grepl("CLON", clonality), 
                Variant_Classification != "Silent")%>%
  group_by(Tumor_Sample_Barcode, clonality)%>%
  dplyr::count()%>%pivot_wider(names_from = clonality, values_from = n)%>%
  mutate(Trial_ID = substr(Tumor_Sample_Barcode, 1, 5),
         ITH = SUBCLONAL/CLONAL)

TCGA_Rx = TCGA_Rx%>%
  dplyr::filter(TREATMENT_TYPE %in% c("Chemotherapy", "Targeted Molecular Therapy"),
                START_DATE>0)

sum(!TCGA_Rx$PATIENT_ID %in% TCGA_clinical$PATIENT_ID)
TCGA_clinical = TCGA_clinical%>%
  mutate(DFS_DAYS = DFS_MONTHS/12*365.25,
         Received_chemo_or_targeted = PATIENT_ID %in% TCGA_Rx$PATIENT_ID,
         Adj_treatment_start = TCGA_Rx$START_DATE[match(PATIENT_ID, TCGA_Rx$PATIENT_ID)],
         total_driver_count = driver_counts_TCGA$total_driver_count[match(PATIENT_ID, driver_counts_TCGA$Tumor_Sample_Barcode)],
         clonal_TMB = ITH_tmb_TCGA$CLONAL[match(PATIENT_ID, ITH_tmb_TCGA$Tumor_Sample_Barcode)],
         subclonal_TMB = ITH_tmb_TCGA$SUBCLONAL[match(PATIENT_ID, ITH_tmb_TCGA$Tumor_Sample_Barcode)],
         ITH = ITH_tmb_TCGA$ITH[match(PATIENT_ID, ITH_tmb_TCGA$Tumor_Sample_Barcode)],
         #Ploidy = seg.mat.copy.list$segments$Ploidy[match(PATIENT_ID, seg.mat.copy.list$segments$SampleID)],
         NStage_binned = substr(PATH_N_STAGE, 1, 2),
         #DFS_DAYS - Adj_treatment_start,
         received_adjuvant_chemo = Received_chemo_or_targeted & Adj_treatment_start <84)%>%
  rowwise()%>%
  mutate(total_TMB = sum(clonal_TMB, subclonal_TMB, na.rm = T),
         total_TMB = if_else(total_TMB == 0, NA, total_TMB))

driver_counts_TCGA$Tumor_Sample_Barcode %in% TCGA_clinical$PATIENT_ID         



TCGA_clinical%>%dplyr::filter(received_adjuvant_chemo == T)%>%pull(Adj_treatment_start)%>%median()


TCGA_clinical$PFS_STATUS%>%table
TCGA_clinical$event = as.numeric(substr(TCGA_clinical$DFS_STATUS, 1, 1))
TCGA_clinical$PFS_MONTHS

TCGA_clinical$STAGE_binned = TCGA_clinical$AJCC_PATHOLOGIC_TUMOR_STAGE%>%sub(".*\\ ", "", .)%>%sub("X|A|B|C", "", .)
TCGA_clinical = TCGA_clinical%>%mutate(STAGE_binned = if_else(STAGE_binned %in% c("I", "II"), "early stage", if_else(STAGE_binned %in% c("III", "IV"), "late stage", NA)))
# general PFS -------------------------------------------------------------
test_df = TCGA_clinical%>%dplyr::filter(grepl("Lum", SUBTYPE))
surv_obj <- Surv(time = test_df$PFS_MONTHS, event = test_df$event)
fit = survfit(surv_obj ~ 1, data = test_df)


survminer::ggsurvplot(fit,
                      data = test_df,
                      risk.table = TRUE,
                      pval = TRUE,
                      palette = nice_cols,
                      title = "Overall disease-free (progression-free) survival"
)


# PFS - TCGA - By Stage ---------------------------------------------------
test_df = TCGA_clinical%>%dplyr::filter(grepl("Lum", SUBTYPE))
surv_obj <- Surv(time = test_df$DFS_MONTHS, event = test_df$event)
fit = survfit(surv_obj ~ STAGE_binned, data = test_df)


survminer::ggsurvplot(fit,
                      data = test_df,
                      risk.table = TRUE,
                      pval = TRUE,
                      palette = nice_cols,
                      title = "DFS by Stage"
)



# STAGE 3 TCGA by Adj Treatment -------------------------------------------
test_df = TCGA_clinical%>%dplyr::filter(STAGE_binned == "late stage",
                                        grepl("Lum", SUBTYPE))
surv_obj <- Surv(time = test_df$DFS_MONTHS, event = test_df$event)
fit = survfit(surv_obj ~ received_adjuvant_chemo, data = test_df)


survminer::ggsurvplot(fit,
                      data = test_df,
                      risk.table = TRUE,
                      pval = TRUE,
                      palette = nice_cols,
                      title = "Stage 3 patients"
)



# CoxPH -------------------------------------------------------------------



TCGA_clinical%>%
  dplyr::filter(PATIENT_ID %in%
  matched_patients_char$PATIENT_ID)%>%
  dplyr::filter(NStage_binned != "NX")%>%
  mutate(metastatic = NStage_binned != "N0" |
           PATH_M_STAGE== "M1" |
          DFS_STATUS == "1:Recurred/Progressed" |
           PFS_STATUS ==  "1:PROGRESSION"
           )%>%
  ggplot(aes(metastatic, subclonal_TMB))+
  geom_violin(aes(fill = metastatic))+
  geom_boxplot(width = .2)+
  geom_jitter(width = 0.02)+
  scale_y_log10()+
  scale_fill_manual(values = c(RS_col, SR_col))+
  theme_classic()+
  stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
  theme_classic(base_size = basesize)+
  theme(legend.position = "none")+
  ylab("Subclonal Mutation Count")

mean(TCGA_clinical$Ploidy [grepl("Lum", TCGA_clinical$SUBTYPE)], na.rm = T)

mean(clinical_data$ploidy, na.rm = T)

  TCGA_clinical%>%
    dplyr::filter(!is.na(STAGE_binned))%>%
    ggplot(aes(STAGE_binned, subclonal_TMB))+
    geom_violin(aes(fill = STAGE_binned))+
    geom_boxplot(width = .2)+
    geom_jitter(width = 0.02)+
    scale_y_log10()+
    scale_fill_manual(values = c(RS_col, SR_col))+
    theme_classic()+
    stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
    theme_classic(base_size = basesize)+
    theme(legend.position = "none")+
    ylab("Subclonal Mutation Count")
  
  

clustering_mat = maf_out%>%
  mutate(mutType = substr(mutType, 3, 5))%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  group_by(Tumor_Sample_Barcode, sample_type, mutType, clonality)%>%
  dplyr::count()%>%
  group_by(Tumor_Sample_Barcode, clonality)%>%
  mutate(sum = sum(n),
         prop = n/sum, 
         mutType = if_else(is.na(mutType), "non_snv", mutType))%>%
  ungroup()%>%
  dplyr::select(prop, mutType, Tumor_Sample_Barcode, clonality)%>%
  pivot_wider(names_from = c(clonality, Tumor_Sample_Barcode), values_from = prop)%>%as.data.frame()
rownames(clustering_mat) = clustering_mat[,1]
clustering_mat = as.matrix(clustering_mat[,2:ncol(clustering_mat)])
clustering_mat[is.na(clustering_mat)] = 0
#clustering_mat = clustering_mat[,colSums(clustering_mat)==1]
clustering_mat = t(clustering_mat)
hc = hclust(dist(clustering_mat), method = "ward.D2")
row_order <- rownames(clustering_mat)[hc$order]
maf_out%>%
  dplyr::filter(sample_type %in% c("FFPE", "ENRICHED"))%>%
  mutate(mutType = substr(mutType, 3, 5),
         mutType = if_else(is.na(mutType), "non_snv", mutType))%>%
  dplyr::filter(grepl("CLON", clonality))%>%
  group_by(Tumor_Sample_Barcode, sample_type, mutType, clonality)%>%
  dplyr::count()%>%
  group_by(Tumor_Sample_Barcode, clonality)%>%
  mutate(sum = sum(n),
         prop = n/sum,
         order = factor(paste(clonality, Tumor_Sample_Barcode, sep = "_"), levels = row_order))%>%
  ggplot(aes(order, prop, fill = mutType))+
  geom_col(width = 1)+
  #scale_y_log10()+
  facet_wrap(~clonality+sample_type, scales = "free_x")+
  theme_classic(base_size = basesize)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top"
  )+
  scale_fill_manual(values = nice_cols)


data_cnv_scores

purity_df = 
  data_cnv_scores%>%
  mutate(SampleIDs = Tumour_ID,
         sample_type = if_else(grepl("ET", SampleIDs), "CKsort", if_else(grepl("_B|FT", SampleIDs), "FFPE", "Other")),
         TrialID = substr(SampleIDs, 1, 5))
FFPE_purity = purity_df%>%dplyr::filter(sample_type == "FFPE")%>%dplyr::filter(Tumour_ID %in% c(maf_out$Tumor_Sample_Barcode[maf_out$SingReg_Prim_v_Met %in% c("PRIMARY", "RECURRENCE_RS")]))  
CK_purity = purity_df%>%dplyr::filter(sample_type == "CKsort", Tumour_ID %in% ck_included_cases_vault)%>%dplyr::select(Purity_CKsort = purity, TrialID, Tumour_ID)  
maf_out$ploidy = purity_df$Ploidy[match(maf_out$Tumor_Sample_Barcode, rownames(purity_df))]

paired_purity = FFPE_purity%>%left_join(CK_purity, by = "TrialID")%>%
  dplyr::filter(!is.na(Purity_CKsort))%>%
  mutate(Purity_FFPE = purity)%>%
  dplyr::select(-purity)%>%
  pivot_longer(cols = c(Purity_FFPE, Purity_CKsort), values_to = "Purity",names_to = "Sample Type")%>%
  mutate(`Sample Type` = if_else(`Sample Type` == "Purity_CKsort", "Purity_RepSeq", `Sample Type`))%>%
  mutate(`Sample Type` = factor(sub("Purity_", "", `Sample Type`), levels = c("FFPE", "RepSeq")))%>%
  ggplot(aes(`Sample Type`, Purity))+
  geom_violin(aes(fill = `Sample Type`), width = 0.8, trim = T, scale = "width")+
  geom_boxplot(width = 0.2, outlier.alpha = 0, varwidth = F)+
  geom_line(aes(`Sample Type`, Purity, group = SampleIDs))+
  geom_point()+
  ylim(0,1)+
  stat_compare_means(paired = T, label.x.npc = 0.4, label = "p.signif", size = label_size, label.y.npc = 0.9)+
  theme_classic(base_size = basesize)+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = c("RepSeq"= RS_col, "FFPE" = SR_col))



# correlation Nuclear Area with Ploidy ------------------------------------

Detections = Detections%>%
  mutate(Sample_ID = substr(Image, 1, 11),
         Trial_ID = substr(Image, 1, 5),
  )

unique(Detections$Sample_ID)

Detections$Sample_ID = sub("HF039_BA006", "HF039_FT003", sub("HF039_BA011", "HF039_FT001", sub("HF044_013.n", "HF044_FT003", sub("HF044_010.n", "HF044_FT001", sub("HF044_009.n", "HF044_FT002", Detections$Sample_ID)))))



data_purity = data_cnv_scores
data_purity_RS = data_purity%>%dplyr::filter(Tumour_ID %in% ck_included_cases_vault)
data_purity = data_purity[grepl("_FT|_B", data_purity$Tumour_ID),]

data_purity$Tumour_ID = sub("_D0.*|_F0.*", "", data_purity$Tumour_ID)





Detections$ploidy_sr = data_purity$ploidy[match(Detections$Sample_ID, data_purity$Tumour_ID)]
Detections$WGD_SR = data_purity$wgd[match(Detections$Sample_ID, data_purity$Tumour_ID)]


Detections$ploidy = data_purity_RS$ploidy[match(Detections$Trial_ID, data_purity_RS$Patient_ID)]
Detections$WGD = data_purity_RS$wgd[match(Detections$Trial_ID, data_purity_RS$Patient_ID)]


plotting_wgd = function(df){
  df+  geom_violin(aes(fill = WGD))+
    geom_boxplot(width = 0.2)+
    geom_jitter(size = pointsize)+
    #geom_text(aes(label = Sample_ID), position = position_jitter())+
    stat_compare_means(label = "p.signif", label.x.npc = .5, label.y.npc = 0.9, size = label_size)+
    theme_classic(base_size = basesize)+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("TRUE" = RS_col, "FALSE" = SR_col))+
    ylab("Nuclear Area (um2)")
}



plotting_ploidy = function(df){
  df+
    geom_point(size = pointsize)+
    # geom_text(aes(label = Sample_ID))+
    geom_smooth(method = "lm", col = "black")+
    stat_cor(size = label_size)+theme_classic(base_size = basesize)+
    ylab("Nuclear Area (um2)")+
    xlab("Ploidy")
}

Detections%>%
  group_by(Sample_ID, WGD, Trial_ID)%>%
  reframe(area = Area.µm.2[Area.µm.2>quantile(Area.µm.2, 0.95, na.rm = T)])%>%
  group_by(WGD)%>%
  reframe(median(area))
pointsize = 2

ploidy_grid = plot_grid(
  Detections%>%
    dplyr::filter(!is.na(WGD_SR))%>%
    group_by(Sample_ID, WGD_SR, Trial_ID)%>%
    reframe(area = Area.µm.2[Area.µm.2>quantile(Area.µm.2, 0.95, na.rm = T)])%>%
    group_by(Sample_ID, WGD_SR, Trial_ID)%>%
    reframe(Nuclear_Area = mean(area))%>%
    mutate(WGD = WGD_SR)%>%
    ggplot(aes(WGD, Nuclear_Area))%>%plotting_wgd()+scale_fill_manual(values = c(SR_col, SR_col))#+
    #ggtitle("Ploidy from FFPE")
  ,
  Detections%>%
    group_by(Sample_ID, ploidy_sr, Trial_ID)%>%
    reframe(area = Area.µm.2[Area.µm.2>quantile(Area.µm.2, 0.95, na.rm = T)])%>%
    group_by(Sample_ID, ploidy_sr, Trial_ID)%>%
    reframe(median_area = median(area))%>%
    mutate(ploidy = ploidy_sr)%>%
    ggplot(aes(ploidy, median_area, col = "black"))%>%plotting_ploidy()+scale_color_manual(values = c(SR_col, SR_col))+
    theme(legend.position = "none"),
  Detections%>%
    group_by(Sample_ID, WGD, Trial_ID)%>%
    reframe(area = Area.µm.2[Area.µm.2>quantile(Area.µm.2, 0.95, na.rm = T)])%>%
    group_by(Sample_ID, WGD, Trial_ID)%>%
    reframe(Nuclear_Area = mean(area))%>%
    ggplot(aes(WGD, Nuclear_Area))%>%plotting_wgd()+scale_fill_manual(values = c(RS_col, RS_col))#+
    #ggtitle("Ploidy from RepSeq")
  ,
  
  Detections%>%
    group_by(Sample_ID, ploidy, Trial_ID)%>%
    reframe(area = Area.µm.2[Area.µm.2>quantile(Area.µm.2, 0.95, na.rm = T)])%>%
    group_by(Sample_ID, ploidy, Trial_ID)%>%
    reframe(median_area = median(area))%>%
    ggplot(aes(ploidy, median_area, col = "black"))%>%plotting_ploidy()+scale_color_manual(values = c(RS_col, RS_col))+
    theme(legend.position = "none"),
  nrow = 1
  
)




plotting_function = function(df, col = "RepSeq"){
  df%>%
    ggplot(aes(maploc, valor, col = col))+
    geom_point(size = 0.1)+
    theme_bw()+
    facet_grid(~chrom, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), strip.background =  element_blank(),
          #strip.text = element_text(face = "bold", size = 12),
          panel.spacing = unit(0, "lines"),
          legend.position = "none"
    )+
    ylab("Variant Allele Log Ratio")
}

facets_snps$chrom[(facets_snps$chrom == 23)]= "X"
facets_snps$chrom = factor(facets_snps$chrom, levels = c(1:22, "X"))

facets_snps_sr$chrom[(facets_snps_sr$chrom == 23)]= "X"
facets_snps_sr$chrom = factor(facets_snps_sr$chrom, levels = c(1:22, "X"))
baf_example_HF270 = plot_grid(facets_snps%>%
            dplyr::filter(!is.na(valor))%>%plotting_function+scale_color_manual(values = RS_col)
          ,
          facets_snps_sr%>%plotting_function+scale_color_manual(values = SR_col),
          nrow = 2)

baf_example_HF270
ggsave(file.path(OUT_DIR, "Figure5h_Example_variantallelelogratio_example.png"))

purity_grid = plot_grid(paired_purity,
purity_across_cohorts,
nrow=1)

plot_grid(purity_grid, ploidy_grid, nrow = 1, rel_widths = c(1:2))
ggsave(file.path(OUT_DIR, "Figure5A_F_SOC_v_RepSeq_purity_ploidy.png"))
