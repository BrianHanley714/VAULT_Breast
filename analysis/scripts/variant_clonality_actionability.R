# Create a plot looking at the proportion of subclonal variants per driver type

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
VARIANTS_MBTCGA = file.path(BASE, "data","variants", "variant_calls_TCGA_MB.txt")
MATCHED_PT_MBTCGA = file.path(BASE, "data","metadata", "matched_patients.txt")
MATCHED_PT_CHAR_MBTCGA = file.path(BASE, "data","metadata", "matched_patients_characteristics.txt")
IDMAP = file.path(BASE, "data","metadata", "tumouridmap_MB.txt")
INCLUDED_PATIENTS = file.path(BASE, "data","metadata", "cases_included.xlsx")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")

# LOAD DATA ---------------------------------------------------------------
rs_patients = read.delim(INCLUDED_PATIENTS)[,1]
vault = read.delim(VARIANTS_VAULT)
mbtcga = read.delim(VARIANTS_MBTCGA)
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
matched_patients = c(matched_patients_char$PATIENT_ID, unique(tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)]))
matched_patients = matched_patients[!is.na(matched_patients)]
clinical_data = read.delim(CLINDATA)


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "custom_filters.R"))
source(file.path(BASE, "src", "rot.lab.R"))
source(file.path(BASE, "src", "annotate_driver_summary.R"))

rs_col = "#54278F"
sr_col = "#33A02Cb7"
sr_col_1 = "#33A02Cb7"
sr_col_2 = "#B2DF8Aff"
gt_col = "grey"
clonal_col = "#33A02Cb7"
subclonal_col= "#54278F"



# COMBINED DATAFRAME ------------------------------------------------------
vault$study = "VAULT"
common_names = Reduce(intersect,list(names(mbtcga), names(vault)))


combined_df = bind_rows(mbtcga%>%dplyr::select(common_names),
                        vault%>%dplyr::select(common_names)
)



# CHISQUARE TEST ----------------------------------------------------------
chisquare = combined_df%>%
  filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%
  mutate(ONCOGENIC = if_else(grepl("Onco", ONCOGENIC_ONCOKB), ONCOGENIC_ONCOKB, "Non_Driver"))%>%
  mutate(Sensitivity_biomarker = if_else(HIGHEST_LEVEL_ONCOKB!= "" &!is.na(HIGHEST_LEVEL_ONCOKB), "Therapeutic_Biomarkers", ONCOGENIC_ONCOKB),
         Sensitivity_biomarker = if_else(grepl("Non_Driver|Therapeutic_Biomarkers", Sensitivity_biomarker), Sensitivity_biomarker, "Driver_Other"),
  )%>%
  filter_rs()%>%
  filter_matched()%>%
  annotate_driver_summary()%>%
  mutate(actionable_summary = if_else(grepl("Biomarker|trial", actionable_summary), "Actionable Driver", if_else(grepl("Other", actionable_summary), "Non-actionable Driver", "Non-Driver")))%>%
  mutate(study = if_else(study %in% c("METABRIC", "TCGA"), "SingReg", "RepSamp"))%>%
  group_by(actionable_summary, study, clonality)%>%
  reframe(count = n())%>%
  group_by(study, actionable_summary)%>%
  mutate(sum = sum(count))%>%
  ungroup()%>%
  pivot_wider(names_from = study, values_from = c(count, sum))%>%
  rowwise()%>%
  mutate(p_value = (chisq.test(matrix(c(count_RepSamp,sum_RepSamp-count_RepSamp, count_SingReg, sum_SingReg-count_SingReg), ncol = 2))$p.value),
         p_label =  ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3))))



# DRAW PLOT  -------------------------------------------------------------
combined_df%>%
  filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%annotate_driver_summary()%>%
  filter_rs()%>%
  filter_matched()%>%
  mutate(study = if_else(study %in% c("METABRIC", "TCGA"), "SingReg", "RepSamp"))%>%
  mutate(actionable_summary = if_else(grepl("Biomarker|trial", actionable_summary), "Actionable Driver", if_else(grepl("Other", actionable_summary), "Non-actionable Driver", "Non-Driver")))%>%
  group_by(actionable_summary, study, clonality)%>%
  reframe(count = n())%>%
  group_by(study, actionable_summary)%>%
  
  mutate(sum = sum(count),
         prop = count/sum,
         p_label = chisquare$p_label[match(actionable_summary, chisquare$actionable_summary)],
         p_label = if_else(study == "SingReg", "", p_label))%>%
  ggplot()+
  geom_col(aes(study, prop, fill = clonality))+
  geom_text(aes(x = study, y = 1.1, label = p_label), nudge_x = 0.5)+
  facet_wrap(~actionable_summary)+
  theme_classic(base_size = 25)+
  scale_fill_manual(values = c("CLONAL" = clonal_col, "SUBCLONAL" = subclonal_col))+
  rot.lab()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 20))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  ylab("Variants (Proportion)")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3E_proportion_variants_clonality_actionability.png"))

