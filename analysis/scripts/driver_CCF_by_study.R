# Compare the cancer cell fraction of drivers between studies

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)

# PATHS -------------------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")
VARIANTS_VAULT = file.path(BASE, "data","variants", "variant_calls_VAULT.txt")
VARIANTS_MBTCGA = file.path(BASE, "data","variants", "variant_calls_TCGA_MB.txt")
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
matched_patients = c(matched_patients_char$PATIENT_ID, tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)])
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


# DRAW PLOT ---------------------------------------------------------------
combined_df%>%
  filter(!is.na(clonality)& clonality != "INDETERMINATE")%>%
  filter_rs()%>%
  filter_matched()%>%annotate_driver_summary()%>%
  mutate(actionable_summary = if_else(grepl("Biomarker|trial", actionable_summary), "Actionable Driver", if_else(grepl("Other", actionable_summary), "Non-actionable Driver", "Non-Driver")))%>%
  mutate(study = if_else(study %in% c("METABRIC", "TCGA"), "SingReg", "RepSamp"))%>%
  filter(actionable_summary == "Actionable Driver")%>%
  filter(clonality == "SUBCLONAL")%>%
  ggplot(aes(study, ccf_expected_copies, fill = study))+
  geom_boxplot()+
  stat_compare_means(size = 6, label.x.npc = 0.25)+
  theme_classic(base_size = 25)+
  ggtitle("Actionable Subclonal Drivers CCF")+
  ylab("CCF")+
  theme(axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_fill_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure3F_actionable_driver_CCF_by_study.png"))

