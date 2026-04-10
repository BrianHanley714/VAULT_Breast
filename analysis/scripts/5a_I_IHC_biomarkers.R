
library(mutSignatures)
library(tidyverse)
library(qs)
library(ggnewscale) 
library(BSgenome.Hsapiens.UCSC.hg19)
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"

dark_grey  = "#4D4D4D"
light_grey = "#D9D9D9"

# PATHS -------------------------------------------------------------------
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
KI67_LOC = file.path(BASE, "data", "image_analysis", "Ki67_pathologist_scores.txt")
FLOW_PATH = file.path(BASE, "data", "fixed_flow", "fixed_flow.txt")
KI67_PATH = file.path(BASE, "data", "image_analysis", "Ki67_pathologist_scores.txt")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")
HF007_LOC = file.path(BASE, "data", "image_analysis", "HF007_image_analysis.txt")


# LOAD DATA ---------------------------------------------------------------
HF007 = read.delim(HF007_LOC)
flow = read.delim(FLOW_PATH)
ki67 = read.delim(KI67_PATH)

clinical_data = read.delim(VAULT_CLIN_PATH)
flow_data = read.delim(FLOW_PATH)
mat_out = read.delim(KI67_LOC)


# LOAD DATA ---------------------------------------------------------------


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))


flow_data = flow_data%>%dplyr::mutate(Trial.ID = substr(Trial_id, 1, 5))
#clinical_data = clinical_data%>%dplyr::mutate(Trial.ID = substr(Trial_id, 1, 5))
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


erprher2 = plot_grid(
  clinical_data%>%
    mutate(ER_low = if_else(ER_low == "TRUE", "Negative", "Positive"))%>%
    ggplot(aes(ER_low, ER_index))+
    geom_violin(aes(fill = ER_low))+
    geom_boxplot(width = 0.2)+
    geom_jitter(width = 0.03)+
    stat_compare_means(label.x.npc = 0.5, label = "p.signif", label.y = 0.9)+
    theme_classic(base_size = basesize)+
    ylab("ER/CK8.18 - FACS")+
    xlab("ER - IHC")+
    scale_fill_manual(values = c(light_grey, light_grey))+
    #scale_y_log10()+
    theme(legend.position = "none",
          axis.title.y = element_text(colour = RS_col),
          axis.title.x = element_text(colour = SR_col)),
  
  clinical_data%>%
    mutate(PR_low = if_else(PR_low == "TRUE", "Negative", "Positive"))%>%
    ggplot(aes(PR_low, PR_index))+
    geom_violin(aes(fill = PR_low))+
    geom_boxplot(width = 0.2)+
    geom_jitter(width = 0.03)+
    stat_compare_means(label.x.npc = 0.5, label = "p.signif", label.y = 0.9)+
    theme_classic(base_size = basesize)+
    ylab("PR/CK8.18 - FACS")+
    xlab("PR - IHC")+
    scale_fill_manual(values = c(light_grey, light_grey))+
    #scale_y_log10()+
    theme(legend.position = "none",
          axis.title.y = element_text(colour = RS_col),
          axis.title.x = element_text(colour = SR_col)),
  
  clinical_data%>%
    ggplot(aes(Her.2_binned, HER2_index))+
    geom_violin(aes(fill = Her.2_binned))+
    geom_boxplot(width = 0.2)+
    geom_jitter(width = 0.03)+
    stat_compare_means(label.x.npc = 0.5, label = "p.signif", label.y = 0.9)+
    theme_classic(base_size = basesize)+
    ylab("HER2/CK8.18 - FACS")+
    xlab("HER2 - IHC")+
    scale_fill_manual(values = c(light_grey, light_grey))+
    #scale_y_log10()+
    theme(legend.position = "none",
          axis.title.y = element_text(colour = RS_col),
          axis.title.x = element_text(colour = SR_col)),
  nrow = 1
)


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))

# GRAPHICAL PARAMETERS ----------------------------------------------------
sr_col = "#33A02Cb7"
rs_col = "#6A3D9Aff"


# DRAW PLOT ---------------------------------------------------------------
# get x-axis order
# mat_out%>%
#   group_by(Trial_ID)%>%
#   reframe(mad = mad(exp_perc))%>%pull(mad)%>%range
# 
# order = 
#   mat_out%>%
#   group_by(filenames)%>%
#   mutate(order = Surface_Area*exp_perc,
#          order = sum(order))%>%
#   arrange(Trial_ID, order)%>%
#   filter(!duplicated(filenames))%>%
#   ungroup()%>%pull(filenames)

ki67_heterogeneity = mat_out%>%  
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         region_mean = sum(order))%>%
  group_by(Trial_ID)%>%
  mutate(max = max(region_mean),
         min = min(region_mean),
         mean = mean(region_mean),
         range = max-min)%>%
  ggplot()+
  geom_jitter(aes(reorder(Trial_ID, mean), exp_perc*100, col = "sub-regional"), alpha = .5, size = 2)+
  geom_point(aes(Trial_ID, region_mean, col = "regional_mean"), alpha = 0.6, size = 3)+
  geom_point(aes(Trial_ID, mean, col = "tumour_mean"), alpha = 0.9, size = 4)+
  scale_color_manual(values = c("tumour_mean" = dark_grey, "regional_mean" = sr_col, "sub-regional" = "grey"))+
  theme_classic(base_size = 20)+
  ylab("Ki67 Positivity (%)")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(colour = SR_col),
        legend.position = c(0.15,0.8))+
  rot.lab()



# GRAPHICAL PARAMETERS ----------------------------------------------------
point_size = 5
base_size = 20
font_size = 8


# ORGANISE DATA -----------------------------------------------------------
flow$Trial_ID = flow$Trial_id
mat_out = left_join(ki67, flow, by = "Trial_ID")
clinical_data = clinical_data%>%select(everything(), Trial_ID = Trial.ID)
flow_out = left_join(flow, clinical_data, by = "Trial_ID")



# DRAW PLOTS --------------------------------------------------------------
plot_1 = 
  mat_out%>%
  group_by(Trial_ID)%>%
  mutate(distinct = n_distinct(filenames))%>%
  group_by(filenames)%>%
  mutate(order = Surface_Area*exp_perc,
         region_mean = sum(order))%>%
  group_by(Trial_ID)%>%
  mutate(max = max(region_mean),
         min = min(region_mean),
         mean = mean(region_mean),
         range = max-min)%>%
  ggplot(aes(mean, KI67.ALL, label = Trial_ID))+
  #geom_text()+
  geom_point(size = point_size)+
  geom_smooth(method = "lm", col = "black")+
  stat_cor(size = font_size)+
  theme_classic(base_size = base_size)+
  ylab("Ki67 (%) Fixed FACS")+
  xlab("Ki67 (%) IHC - Tumour Mean")+
  theme(legend.position = "none",
        axis.title.y = element_text(colour = RS_col),
        axis.title.x = element_text(colour = SR_col))



plot_2 =
  flow_out%>%
  dplyr::filter(!is.na(Pathology.Grading.Components_M))%>%
  ggplot(aes(as.character(Pathology.Grading.Components_M), KI67.ALL.x))+
  #geom_boxplot(outlier.alpha = 0)+
  #geom_jitter(size = point_size)+
  geom_violin(fill = "lightgrey", width= 1.2)+
  geom_boxplot(width = 0.5, outlier.alpha = 0)+
  geom_jitter(size = point_size/2, width = 0.11)+
  theme_classic(base_size = base_size)+
  xlab("Mitotic Rate - H&E")+
  ylab("Ki67 (%) Fixed FACS")+
  stat_compare_means(comparisons = list(c("1", "2"),
                                        c("1", "3"),
                                        c("3", "2")),label = "p.signif", size = font_size)+
  scale_y_log10()+
  theme(legend.position = "none",
        axis.title.y = element_text(colour = RS_col),
        axis.title.x = element_text(colour = SR_col))


plot_3 = 
  flow_out%>%
  mutate(Recurrence.score_GEP = if_else(Recurrence.score_GEP == 0, NA, Recurrence.score_GEP))%>%
  mutate(label = if_else(Trial_ID == "HF007", "HF007", ""))%>%
  dplyr::filter(!is.na(Recurrence.score_GEP))%>%
  ggplot(aes(Recurrence.score_GEP, KI67.ALL.y, label = label))+
  geom_point(size = point_size, aes(col = GEP_performed))+
  geom_smooth(method = "lm", col = "black", se = F)+
  theme_classic(base_size = base_size)+
  geom_text(nudge_x = 5, fontface = "bold", size = font_size*0.75)+
  stat_cor(size = font_size)+
  scale_color_manual(values = c(light_grey, dark_grey))+
  theme(legend.title = element_blank(),
        legend.position = c(0.3,0.6),
        axis.title.y = element_text(colour = RS_col),
        axis.title.x = element_text(colour = SR_col),
        legend.text = element_text(size = font_size*2))+
  ylab("Ki67 (%) Fixed FACS")+
  xlab("Recurrence Score")

plot_grid(
erprher2,
plot_grid(plot_1, plot_2, plot_3, ncol = 3, label_size = 25),
nrow = 2
)
ggsave(file.path(OUT_DIR, "Figure7_ER_PR_HER2_KI67_IHC_FLOW_Correlations.png"))

ki67_heterogeneity
ggsave(file.path(OUT_DIR, "Figure7_Ki67_heterogeneity.png"))
HF007%>%
  
  dplyr::filter(Name == "Tumor")%>%
  mutate(OncotypeDX_performed = if_else(block == "BA015", TRUE, FALSE),
         slice = case_when(block == "BA010"~"Slice_7",
                           block == "BA012"~"Slice_8",
                           block == "BA014"~"Slice_10",
                           block == "BA015"~"Slice_11",
                           block == "BA016"~"Slice_12",
                           block == "BA018"~"Slice_13",
         ),
         order = as.numeric(sub("BA0", "", block)))%>%
  ggplot(aes(reorder(slice, order), Tumour_cells..Positive..))+
  geom_boxplot(outlier.alpha = 0, aes(fill = OncotypeDX_performed))+
  geom_jitter(size = 3)+
  theme_classic(base_size = 20)+
  stat_compare_means(label.x = 3, label.y = 40, label = "p.signif", size = 10)+
  scale_fill_manual(values = c("TRUE" = light_grey, "FALSE" = dark_grey))+
  ylab("Ki67+ Cancer Cells (%)")+
  theme(axis.title.x = element_blank(), 
        legend.position = "top")

ggsave(file.path(OUT_DIR, "Figure7_HF007_heterogeneity.png"))
