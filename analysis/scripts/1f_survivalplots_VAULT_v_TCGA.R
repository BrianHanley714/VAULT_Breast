
rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(survival)
library(survminer)
library(MatchIt)
library(ggnewscale)
library(openxlsx)
library(qs)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
TCGA_CLIN_PATH = file.path(BASE, "data", "metadata", "whole_TCGA_Breast_metadata.txt")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")
VAULT_CLIN_PATH = file.path(BASE, "data", "metadata", "clinical_data.txt")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")
# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))



# LOAD DATA ---------------------------------------------------------------
survival_tcga = read.delim(TCGA_CLIN_PATH)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]
survival_vault = read.delim(VAULT_CLIN_PATH)
maf_out = qread(VAULT_VARIANTS_PATH)
# options -----------------------------------------------------------------
basesize = 20
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"









# ORGANISE DATA -----------------------------------------------------------

survival_vault$ICD10 = case_when(survival_vault$Histological.Subtype == "IDC_NST" ~ "8500/3",
                                 survival_vault$Histological.Subtype == "ILC"~ "8520/3" ,
                                 survival_vault$Histological.Subtype == "MixedILC/IDC"~"8522/3",
                                 survival_vault$Histological.Subtype == "MucinousCarcinoma"~"8480/3")

survival_vault$subtype = survival_vault$Receptor.Subtype
survival_vault$subtype = case_when(survival_vault$Receptor.Subtype == "ER+/HER2-" ~"Luminal",
                          survival_vault$Receptor.Subtype == "HER2+" ~"HER2-overexpressing",
                          survival_vault$Receptor.Subtype == "TN" ~"TN")


survival_tcga = survival_tcga%>%
  dplyr::filter(!is.na(molsubtype)& !is.na(ICD0))%>%
  dplyr::filter(ICD0 %in% unique(survival_vault$ICD10))






survival_out = bind_rows(
  survival_tcga%>%
    mutate(status = if_else(DFS_STATUS == "0:DiseaseFree", 0, 1),
           daystolastfollow = DFS_MONTHS*30.4375,
           TCGA = "TCGA")%>%
    dplyr::filter(!is.na(daystolastfollow))%>%
    dplyr::select(PATIENT_ID,  status, daystolastfollow, study = "TCGA", ICD10 = ICD0, molsubtype, Age = AGE),
  
  survival_vault%>%
    dplyr::filter(Trial.ID %in% c(substr(ck_included_cases_vault, 1, 5)),recurrence_current_sample_is_recurence == "No")%>%
    mutate(status = if_else(is.na(days_to_first_recurrence), 0, 1),
           daystolastfollow = if_else(is.na(days_to_first_recurrence), daystolastfollow, days_to_first_recurrence),
           VAULT = "VAULT")%>%
    dplyr::select(PATIENT_ID= Trial.ID, status, daystolastfollow, study = "VAULT", ICD10, molsubtype = subtype, Age)
)



df = as.data.frame(survival_out%>%dplyr::filter(study != "METABRIC", !is.na(ICD10), !is.na(molsubtype))%>%mutate(study = if_else(study == "VAULT", 1, 0)))


m.out0 <- matchit(study ~ ICD10 + molsubtype+Age,
                  data = df,
                  method = NULL,
                  distance = "glm")

out = matchit(study ~ ICD10 + molsubtype,
              data = df,
              method = "optimal",replace = F, ratio = 5,
              distance = "glm")


out = match.data(out)

out = out%>%
  mutate(study = if_else(study == 1, "VAULT", "TCGA"))
fit = survfit(Surv(daystolastfollow, status) ~ study, data = out)



ggsurvplot(
  fit,
  data = out,
  risk.table = F,       
  pval = TRUE,             
  conf.int = F,          
  legend.title = "Disease Free Survival",
  palette = c(SR_col, RS_col, "red"),
  ggtheme = theme_classic(base_size = basesize) 
)


# Write Figure 1F ---------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure1F_survival_TCGA_V_VAULT.png"))


# Organise data for Supplementary Figure ----------------------------------

# Function
shade_between_white <- function(col, n = 5) {
  target = grDevices::col2rgb(col) / 255
  
  white = c(1, 1, 1)
  
  shades = t(sapply(seq(0, 1, length.out = n), function(a) {
    (1 - a) * white + a * target
  }))
  
  grDevices::rgb(shades[,1], shades[,2], shades[,3])
}

blue_shades = shade_between_white("navy", 4)
red_shades = shade_between_white("maroon", 4)

cases_w_biomarkers = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in%ck_included_cases_vault,
                grepl("Biomarker", actionable_summary))%>%
  pull(Trial_ID)%>%unique

cases_w_trial_inclusion = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in%ck_included_cases_vault,
                grepl("trial_inclusion_criterion", actionable_summary))%>%
  pull(Trial_ID)%>%unique

gap = 40
recurrence.label.x = 1000
legend.label.y=0.3
survival_vault
textsize = 3
square_size = 2


# Output the supplementary figure -------------------------------------------------------


survival_vault%>%
  dplyr::filter(Trial.ID %in% c(substr(ck_included_cases_vault, 1, 5)))%>%
  mutate(status = if_else(is.na(days_to_first_recurrence), 0, 1),
         daystolastfollow = if_else(is.na(days_to_first_recurrence), daystolastfollow, days_to_first_recurrence),
         VAULT = "VAULT")%>%
  dplyr::select(everything(), PATIENT_ID= Trial.ID, study = "VAULT", molsubtype = subtype)%>%
  #dplyr::filter(study == "VAULT")%>%
  #dplyr::filter(PATIENT_ID %in% c(substr(ck_included_cases_vault, 1, 5)))%>%
  mutate(Trial.ID = PATIENT_ID)%>%
  #left_join(survival_vault[,!names(survival_vault) %in% names(survival_out)], by = "Trial.ID")%>%
  mutate(status = if_else(status== 1, "X", NA), 
         T.Stage = substr(T.Stage, 1, 2),
         N.stage = substr(N.stage, 1, 2),
         `Biomarker Identified` = if_else(Trial.ID %in% cases_w_biomarkers, "Yes", "No"),
         `Trial Inclusion Identified` = if_else(Trial.ID %in% cases_w_trial_inclusion, "Yes", "No")
  )%>%
  ggplot(aes(daystolastfollow, reorder(PATIENT_ID, daystolastfollow)))+
  geom_col()+
  geom_text(aes(label = status),col = "red", fontface = "bold", size = text_size/2)+
  geom_point(aes(-gap, reorder(PATIENT_ID, daystolastfollow), fill = T.Stage), shape = 22, size = square_size)+
  scale_fill_manual(values= c("T1" = blue_shades[1], "T2" = blue_shades[2], "T3" = blue_shades[3], "T4" = blue_shades[4]))+
  new_scale_fill()+
  geom_point(aes(-gap*2, reorder(PATIENT_ID, daystolastfollow), fill = N.stage), shape = 22, size = square_size,inherit.aes = FALSE)+
  scale_fill_manual(values= c("N0" = red_shades[1], "N1" = red_shades[2], "N2" = red_shades[3], "N3" = red_shades[4]))+
  new_scale_fill()+
  geom_point(aes(-gap*3, reorder(PATIENT_ID, daystolastfollow), fill = `Biomarker Identified`), shape = 22, size = square_size)+
  geom_point(aes(-gap*4, reorder(PATIENT_ID, daystolastfollow), fill = `Trial Inclusion Identified`), shape = 22, size = square_size)+
  scale_fill_manual(values= c("Yes" = RS_col, "No" = "white"))+
  theme_classic(base_size = basesize*0.75)+
  annotate("text", x = -gap, y = length(ck_included_cases_vault)+1, label = "T Stage", angle = 90, size = textsize, hjust = 0)+
  annotate("text", x = -gap*2, y = length(ck_included_cases_vault)+1, label = "N Stage", angle = 90, size = textsize, hjust = 0)+
  annotate("text", x = -gap*3, y = length(ck_included_cases_vault)+1, label = "Biomarker", angle = 90, size = textsize, hjust = 0)+
  annotate("text", x = -gap*4, y = length(ck_included_cases_vault)+1, label = "Trial IC", angle = 90, size = textsize, hjust = 0)+
  annotate("text", x = recurrence.label.x, y = length(ck_included_cases_vault)+3, label = "X", fontface = "bold",col = "red", size = text_size/2)+
  annotate("text", x = recurrence.label.x+50, y = length(ck_included_cases_vault)+3, label = "recurrence",col = "black", hjust = 0, size = text_size/2)+
  geom_point(aes(x = 1, y = length(ck_included_cases_vault)+12), col= "white")+
  xlab("Days to first recurrence")+
  #  guides(fill = guide_legend(nrow = 1))+
  theme(axis.title.y = element_blank(),legend.title = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y  = element_text(margin = margin(r = -15), size = textsize*2.5),
        axis.ticks.y = element_blank(),
        legend.position = c(0.85, legend.label.y))


# Write supplementary figure ----------------------------------------------

ggsave(file.path(OUT_DIR, "SupplementaryFigure_Disease_free_survival_features.png"))

